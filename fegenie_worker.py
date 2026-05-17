#!/usr/bin/env python3
"""
FeGenie S3 worker for the Amplify React app.

What it does
------------
1. Scans an input S3 bucket for uploaded FeGenie jobs:
     s3://<input-bucket>/fegenie-<slug>/manifest-<slug>.txt
2. Downloads the job folder to a local work directory.
3. Validates that uploaded genomes are either FASTA contigs OR GenBank files.
4. Normalizes input file extensions into a clean run directory.
5. Runs FeGenie.
6. Writes frontend-ready results to:
     s3://<results-bucket>/fegenie-<slug>/

Frontend-required result names
------------------------------
The React results viewer fetches exactly:
  - FeGenie-geneSummary.csv
  - FeGenie-heatmap-data.csv
  - fegenie-report.html

This script makes sure those names exist. It also uploads:
  - status.json
  - run.log
  - raw-results.tar.gz

Typical cron usage
------------------
Run every few minutes on the FeGenie worker host:

  /usr/bin/python3 /opt/fegenie/fegenie_worker.py \
    --input-bucket midauthorbio-fegenie-input \
    --results-bucket midauthorbio-fegenie-results \
    --work-root /data/fegenie-worker \
    --fegenie-bin /home/ark/MAB/bin/FeGenie/FeGenie.py \
    --command-prefix "conda run -n fegenie"

If FeGenie.py is already on PATH in the active environment, omit --command-prefix
and set --fegenie-bin FeGenie.py.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from typing import Iterable

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError


FASTA_EXTS = {".fa", ".fasta", ".fna"}
GENBANK_EXTS = {".gb", ".gbk", ".gbff", ".gbf"}
SKIP_INPUT_NAMES = {"form-data.txt"}
MANIFEST_RE = re.compile(r"^manifest-(?P<slug>[A-Za-z0-9]+)\.txt$")
PREFIX_RE = re.compile(r"^fegenie-(?P<slug>[A-Za-z0-9]+)/$")


@dataclass
class JobManifest:
    slug: str
    job_name: str = ""
    submitter_name: str = ""
    submitter_email: str = ""
    analysis_mode: str = "single_genome"
    threads: int = 4
    out: str = "fegenie_out"
    cluster_distance: int | None = None
    inflation: str | None = None
    raw_options: dict[str, str] = field(default_factory=dict)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def configure_logging(log_path: Path, verbose: bool = False) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [
        logging.FileHandler(log_path),
        logging.StreamHandler(sys.stdout),
    ]
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers,
    )


def make_s3(region: str):
    # Important: use the regional endpoint, otherwise browser presigned PUTs may
    # CORS-fail after an S3 redirect from s3.amazonaws.com to the bucket region.
    return boto3.client(
        "s3",
        region_name=region,
        endpoint_url=f"https://s3.{region}.amazonaws.com",
        config=Config(signature_version="s3v4", s3={"addressing_style": "virtual"}),
    )


def make_ses(region: str):
    return boto3.client("ses", region_name=region)


def s3_key_exists(s3, bucket: str, key: str) -> bool:
    try:
        s3.head_object(Bucket=bucket, Key=key)
        return True
    except ClientError as e:
        code = e.response.get("Error", {}).get("Code")
        if code in {"404", "NoSuchKey", "NotFound"}:
            return False
        raise


def put_json(s3, bucket: str, key: str, payload: dict) -> None:
    s3.put_object(
        Bucket=bucket,
        Key=key,
        Body=json.dumps(payload, indent=2, sort_keys=True).encode("utf-8"),
        ContentType="application/json",
    )


def upload_file(s3, bucket: str, key: str, path: Path, content_type: str | None = None) -> None:
    extra = {}
    if content_type:
        extra["ContentType"] = content_type
    s3.upload_file(str(path), bucket, key, ExtraArgs=extra)


def list_job_prefixes(s3, input_bucket: str) -> list[str]:
    paginator = s3.get_paginator("list_objects_v2")
    prefixes: list[str] = []
    for page in paginator.paginate(Bucket=input_bucket, Prefix="fegenie-", Delimiter="/"):
        for item in page.get("CommonPrefixes", []):
            prefix = item["Prefix"]
            if PREFIX_RE.match(prefix):
                prefixes.append(prefix)
    return sorted(prefixes)


def find_manifest_key(s3, input_bucket: str, prefix: str) -> str | None:
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=input_bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            name = Path(obj["Key"]).name
            if MANIFEST_RE.match(name):
                return obj["Key"]
    return None


def download_prefix(s3, bucket: str, prefix: str, dest: Path) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    downloaded: list[Path] = []
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith("/"):
                continue
            rel = Path(key).relative_to(prefix)
            local = dest / rel
            local.parent.mkdir(parents=True, exist_ok=True)
            s3.download_file(bucket, key, str(local))
            downloaded.append(local)
            logging.info("Downloaded s3://%s/%s -> %s", bucket, key, local)
    return downloaded


def parse_manifest(path: Path, fallback_slug: str) -> JobManifest:
    manifest = JobManifest(slug=fallback_slug)
    in_options = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            raw = line.rstrip("\n")
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped == "# Options":
                in_options = True
                continue
            if stripped.startswith("#"):
                if stripped != "# Options":
                    in_options = False
                continue

            if ":" not in stripped:
                continue

            key, value = stripped.split(":", 1)
            key = key.strip().lstrip("-")
            value = value.strip().strip('"')

            if key == "job_slug":
                manifest.slug = value or fallback_slug
            elif key == "job_name":
                manifest.job_name = value
            elif key == "submitter_name":
                manifest.submitter_name = value
            elif key == "submitter_email":
                manifest.submitter_email = value
            elif in_options:
                manifest.raw_options[key] = value

    # Current frontend option keys.
    mode = manifest.raw_options.get("analysis_mode")
    if mode:
        manifest.analysis_mode = mode
    if manifest.raw_options.get("t"):
        manifest.threads = int(float(manifest.raw_options["t"]))
    if manifest.raw_options.get("out"):
        manifest.out = manifest.raw_options["out"]
    if manifest.raw_options.get("d"):
        manifest.cluster_distance = int(float(manifest.raw_options["d"]))
    if manifest.raw_options.get("inflation"):
        manifest.inflation = manifest.raw_options["inflation"]

    return manifest


def safe_output_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def prepare_bins(job_dir: Path, bins_dir: Path) -> tuple[str, list[Path]]:
    """Normalize uploaded input files for FeGenie.

    FeGenie wants one -bin_ext. The web app lets users upload common FASTA or
    GenBank extensions, so we copy files to a clean directory and normalize
    extensions to either .fa or .gbk.
    """
    bins_dir.mkdir(parents=True, exist_ok=True)
    inputs: list[Path] = []

    for path in sorted(job_dir.iterdir()):
        if not path.is_file():
            continue
        lower = path.name.lower()
        if lower in SKIP_INPUT_NAMES or lower.startswith("manifest-"):
            continue
        suffix = path.suffix.lower()
        if suffix in FASTA_EXTS or suffix in GENBANK_EXTS:
            inputs.append(path)

    if not inputs:
        raise RuntimeError("No FASTA or GenBank input files found in the job folder.")

    suffixes = {p.suffix.lower() for p in inputs}
    has_fasta = any(s in FASTA_EXTS for s in suffixes)
    has_genbank = any(s in GENBANK_EXTS for s in suffixes)
    if has_fasta and has_genbank:
        raise RuntimeError("Mixed FASTA and GenBank uploads are not supported in one FeGenie job.")

    target_ext = ".gbk" if has_genbank else ".fa"
    bin_ext = target_ext.lstrip(".")

    seen: set[str] = set()
    normalized: list[Path] = []
    for src in inputs:
        base = safe_output_name(src.stem)
        name = f"{base}{target_ext}"
        i = 2
        while name in seen:
            name = f"{base}_{i}{target_ext}"
            i += 1
        seen.add(name)
        dest = bins_dir / name
        shutil.copy2(src, dest)
        normalized.append(dest)

    return bin_ext, normalized


def build_fegenie_command(args, manifest: JobManifest, bins_dir: Path, out_dir: Path, bin_ext: str) -> list[str]:
    cmd: list[str] = []
    if args.command_prefix:
        cmd.extend(args.command_prefix.split())
    cmd.extend([
        args.fegenie_bin,
        "-bin_dir",
        str(bins_dir),
        "-bin_ext",
        bin_ext,
        "-out",
        str(out_dir),
        "-t",
        str(manifest.threads),
    ])

    if bin_ext == "gbk":
        cmd.append("--gbk")
    if manifest.analysis_mode == "metagenomic":
        cmd.append("--meta")
    if manifest.cluster_distance is not None:
        cmd.extend(["-d", str(manifest.cluster_distance)])
    if manifest.inflation:
        cmd.extend(["-inflation", str(manifest.inflation)])

    return cmd


def run_command(cmd: list[str], log_handle) -> None:
    logging.info("Running command: %s", " ".join(cmd))
    proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}: {' '.join(cmd)}")


def first_existing(out_dir: Path, names: Iterable[str]) -> Path | None:
    for name in names:
        p = out_dir / name
        if p.exists():
            return p
    return None


def normalize_gene_summary(src: Path, dest: Path) -> None:
    """Remove duplicate header rows and ensure frontend-required columns are present."""
    with src.open("r", newline="", encoding="utf-8", errors="replace") as inp:
        rows = list(csv.reader(inp))
    if not rows:
        raise RuntimeError(f"{src} is empty.")

    # FeGenie variants sometimes include a shorter header then a fuller header.
    # Prefer the header with the most expected columns.
    expected = [
        "category",
        "genome/assembly",
        "orf",
        "HMM",
        "bitscore",
        "bitscore_cutoff",
        "clusterID",
        "heme_c_binding_motifs",
        "heme_b_binding_motifs",
        "hematite_binding_motifs",
        "protein_sequence",
    ]

    header_idx = 0
    best_score = -1
    for i, row in enumerate(rows[:5]):
        score = sum(1 for col in expected if col in row)
        if row and row[0] == "category" and score > best_score:
            header_idx = i
            best_score = score

    header = rows[header_idx]
    data = rows[header_idx + 1 :]
    cleaned = [r for r in data if r and r[0] != "category"]

    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(header)
        writer.writerows(cleaned)


def normalize_heatmap(src: Path, dest: Path) -> None:
    with src.open("r", newline="", encoding="utf-8", errors="replace") as inp:
        reader = csv.reader(inp)
        rows = [r for r in reader if r]
    if not rows:
        raise RuntimeError(f"{src} is empty.")
    if rows[0][0] != "X":
        rows[0][0] = "X"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerows(rows)


def write_fallback_report(dest: Path, manifest: JobManifest, heatmap_csv: Path, gene_csv: Path) -> None:
    dest.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>FeGenie Report {manifest.slug}</title>
  <style>
    body {{ font-family: Inter, system-ui, sans-serif; margin: 2rem; line-height: 1.5; color: #241f1c; }}
    code {{ background: #f4eee8; padding: 0.15rem 0.35rem; border-radius: 4px; }}
    .card {{ border: 1px solid #dfd4ca; border-radius: 12px; padding: 1rem; max-width: 760px; }}
  </style>
</head>
<body>
  <h1>FeGenie report</h1>
  <div class="card">
    <p><strong>Job:</strong> <code>{manifest.slug}</code></p>
    <p><strong>Job name:</strong> {manifest.job_name or "unnamed"}</p>
    <p><strong>Analysis mode:</strong> {manifest.analysis_mode}</p>
    <p>Interactive visualizations are available in the FeGenie Web Studio results page.</p>
    <ul>
      <li>{heatmap_csv.name}</li>
      <li>{gene_csv.name}</li>
    </ul>
  </div>
</body>
</html>
""",
        encoding="utf-8",
    )


def maybe_generate_report(args, out_dir: Path, report_dest: Path, manifest: JobManifest, heatmap_csv: Path, gene_csv: Path, log_handle) -> None:
    existing = first_existing(out_dir, ["fegenie-report.html", "FeGenie-report.html", "FeGenie_Report.html"])
    if existing:
        shutil.copy2(existing, report_dest)
        return

    if args.report_script:
        script = Path(args.report_script)
        if script.exists():
            cmd: list[str] = []
            if args.command_prefix:
                cmd.extend(args.command_prefix.split())
            cmd.extend([str(script), "-o", str(report_dest), str(heatmap_csv)])
            run_command(cmd, log_handle)
            if report_dest.exists():
                return

    write_fallback_report(report_dest, manifest, heatmap_csv, gene_csv)


def make_tarball(source_dir: Path, tar_gz: Path) -> None:
    with tarfile.open(tar_gz, "w:gz") as tar:
        tar.add(source_dir, arcname=source_dir.name)


def send_ses_email(args, manifest: JobManifest, subject: str, body: str) -> None:
    if not args.ses_sender or not manifest.submitter_email:
        return

    ses = make_ses(args.region)
    msg = MIMEMultipart()
    msg["From"] = args.ses_sender
    msg["To"] = manifest.submitter_email
    msg["Subject"] = subject
    if args.ses_bcc:
        msg["Bcc"] = args.ses_bcc
    msg.attach(MIMEText(body, "plain"))

    destinations = [manifest.submitter_email]
    if args.ses_bcc:
        destinations.append(args.ses_bcc)

    ses.send_raw_email(
        Source=args.ses_sender,
        Destinations=destinations,
        RawMessage={"Data": msg.as_string()},
    )


def process_job(args, s3, prefix: str) -> None:
    slug_match = PREFIX_RE.match(prefix)
    if not slug_match:
        return
    slug = slug_match.group("slug")
    result_prefix = f"fegenie-{slug}/"
    status_key = f"{result_prefix}status.json"

    if s3_key_exists(s3, args.results_bucket, status_key) and not args.force:
        logging.info("Skipping %s; result status already exists.", prefix)
        return

    manifest_key = find_manifest_key(s3, args.input_bucket, prefix)
    if not manifest_key:
        logging.info("Skipping %s; no manifest file found yet.", prefix)
        return

    job_root = Path(args.work_root) / f"fegenie-{slug}"
    if job_root.exists():
        shutil.rmtree(job_root)
    job_dir = job_root / "input"
    bins_dir = job_root / "bins"
    out_dir = job_root / "fegenie_out"
    final_dir = job_root / "frontend_results"
    log_path = job_root / "run.log"
    job_root.mkdir(parents=True, exist_ok=True)

    put_json(
        s3,
        args.results_bucket,
        status_key,
        {"slug": slug, "state": "running", "started_at": utc_now(), "input_prefix": prefix},
    )

    try:
        download_prefix(s3, args.input_bucket, prefix, job_dir)
        manifest_path = job_dir / Path(manifest_key).name
        manifest = parse_manifest(manifest_path, fallback_slug=slug)
        bin_ext, normalized_inputs = prepare_bins(job_dir, bins_dir)

        with log_path.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"FeGenie worker started {utc_now()}\n")
            log_handle.write(f"Slug: {slug}\n")
            log_handle.write(f"Input files: {[p.name for p in normalized_inputs]}\n")
            cmd = build_fegenie_command(args, manifest, bins_dir, out_dir, bin_ext)
            run_command(cmd, log_handle)

            final_dir.mkdir(parents=True, exist_ok=True)
            heatmap_src = first_existing(out_dir, ["FeGenie-heatmap-data.csv", "heatmap-data.csv"])
            gene_src = first_existing(out_dir, ["FeGenie-geneSummary.csv", "geneSummary.csv"])
            if not heatmap_src:
                raise RuntimeError("Could not find FeGenie-heatmap-data.csv in output directory.")
            if not gene_src:
                raise RuntimeError("Could not find FeGenie-geneSummary.csv in output directory.")

            heatmap_dest = final_dir / "FeGenie-heatmap-data.csv"
            gene_dest = final_dir / "FeGenie-geneSummary.csv"
            report_dest = final_dir / "fegenie-report.html"
            normalize_heatmap(heatmap_src, heatmap_dest)
            normalize_gene_summary(gene_src, gene_dest)
            maybe_generate_report(args, out_dir, report_dest, manifest, heatmap_dest, gene_dest, log_handle)

        tar_path = job_root / "raw-results.tar.gz"
        make_tarball(out_dir, tar_path)

        upload_file(s3, args.results_bucket, f"{result_prefix}FeGenie-heatmap-data.csv", final_dir / "FeGenie-heatmap-data.csv", "text/csv")
        upload_file(s3, args.results_bucket, f"{result_prefix}FeGenie-geneSummary.csv", final_dir / "FeGenie-geneSummary.csv", "text/csv")
        upload_file(s3, args.results_bucket, f"{result_prefix}fegenie-report.html", final_dir / "fegenie-report.html", "text/html")
        upload_file(s3, args.results_bucket, f"{result_prefix}raw-results.tar.gz", tar_path, "application/gzip")
        upload_file(s3, args.results_bucket, f"{result_prefix}run.log", log_path, "text/plain")

        result_url = f"{args.app_url.rstrip('/')}/#fegenie/results-{slug}" if args.app_url else ""
        put_json(
            s3,
            args.results_bucket,
            status_key,
            {
                "slug": slug,
                "state": "complete",
                "completed_at": utc_now(),
                "result_prefix": result_prefix,
                "result_url": result_url,
                "files": [
                    "FeGenie-geneSummary.csv",
                    "FeGenie-heatmap-data.csv",
                    "fegenie-report.html",
                    "raw-results.tar.gz",
                    "run.log",
                ],
            },
        )

        if args.ses_sender and manifest.submitter_email:
            body = (
                f"Hi {manifest.submitter_name or 'there'},\n\n"
                f"Your FeGenie job is complete.\n\n"
                f"Result code: {slug}\n"
                f"Results URL: {result_url or '(configure --app-url to include this)'}\n\n"
                "Results are retained according to the S3 lifecycle rule on the results bucket.\n"
            )
            send_ses_email(args, manifest, "Your FeGenie results are ready", body)

        logging.info("Completed job %s -> s3://%s/%s", slug, args.results_bucket, result_prefix)

    except Exception as e:
        logging.exception("Job %s failed", slug)
        put_json(
            s3,
            args.results_bucket,
            status_key,
            {"slug": slug, "state": "failed", "failed_at": utc_now(), "error": str(e)},
        )
        if log_path.exists():
            upload_file(s3, args.results_bucket, f"{result_prefix}run.log", log_path, "text/plain")
        try:
            manifest = parse_manifest(job_dir / Path(manifest_key).name, fallback_slug=slug)
            send_ses_email(
                args,
                manifest,
                "FeGenie job failed",
                f"FeGenie job {slug} failed.\n\nError:\n{e}\n\nCheck s3://{args.results_bucket}/{result_prefix}run.log",
            )
        except Exception:
            logging.exception("Unable to send failure notification.")
        if not args.continue_on_error:
            raise
    finally:
        if args.clean and job_root.exists():
            shutil.rmtree(job_root)


def acquire_lock(lock_path: Path):
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
    except FileExistsError:
        raise RuntimeError(f"Lock already exists: {lock_path}")


def release_lock(lock_path: Path):
    try:
        lock_path.unlink()
    except FileNotFoundError:
        pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Poll S3 for FeGenie jobs, run FeGenie, and publish frontend-ready results.")
    parser.add_argument("--input-bucket", default=os.getenv("FEGENIE_INPUT_BUCKET", "midauthorbio-fegenie-input"))
    parser.add_argument("--results-bucket", default=os.getenv("FEGENIE_RESULTS_BUCKET", "midauthorbio-fegenie-results"))
    parser.add_argument("--region", default=os.getenv("AWS_REGION", "us-east-2"))
    parser.add_argument("--work-root", default=os.getenv("FEGENIE_WORK_ROOT", "/tmp/fegenie-worker"))
    parser.add_argument("--fegenie-bin", default=os.getenv("FEGENIE_BIN", "FeGenie.py"))
    parser.add_argument("--command-prefix", default=os.getenv("FEGENIE_COMMAND_PREFIX", ""), help='Optional prefix, e.g. "conda run -n fegenie"')
    parser.add_argument("--report-script", default=os.getenv("FEGENIE_REPORT_SCRIPT", ""), help="Optional legacy fegenie_report.py path.")
    parser.add_argument("--app-url", default=os.getenv("FEGENIE_APP_URL", ""), help="Amplify app base URL for optional notifications.")
    parser.add_argument("--ses-sender", default=os.getenv("FEGENIE_SES_SENDER", ""), help="Verified SES sender. If omitted, no email is sent.")
    parser.add_argument("--ses-bcc", default=os.getenv("FEGENIE_SES_BCC", ""))
    parser.add_argument("--once", action="store_true", help="Run one scan and exit.")
    parser.add_argument("--interval", type=int, default=300, help="Polling interval seconds when not using --once.")
    parser.add_argument("--force", action="store_true", help="Reprocess jobs even if result status.json exists.")
    parser.add_argument("--clean", action="store_true", help="Delete local work directory after each job.")
    parser.add_argument("--continue-on-error", action="store_true", help="Continue polling other jobs after a failure.")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--lock-file", default=os.getenv("FEGENIE_LOCK_FILE", "/tmp/fegenie-worker.lock"))
    parser.add_argument("--log-file", default=os.getenv("FEGENIE_WORKER_LOG", "/tmp/fegenie-worker.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(Path(args.log_file), verbose=args.verbose)
    lock_path = Path(args.lock_file)

    try:
        acquire_lock(lock_path)
    except RuntimeError as e:
        logging.warning("%s", e)
        return 0

    s3 = make_s3(args.region)
    try:
        while True:
            prefixes = list_job_prefixes(s3, args.input_bucket)
            logging.info("Found %d FeGenie job prefix(es).", len(prefixes))
            for prefix in prefixes:
                process_job(args, s3, prefix)
            if args.once:
                break
            time.sleep(args.interval)
    finally:
        release_lock(lock_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
