import csv
from pathlib import Path
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import fegenie_runtime as runtime
import fegenie_worker as worker


def feature_record(dna='ATGAAATAA', strand=1, **qualifiers):
    record = SeqRecord(Seq(dna), id='contig.1', name='contig')
    record.annotations['molecule_type'] = 'DNA'
    feature = SeqFeature(SimpleLocation(0, len(dna), strand=strand), type='CDS',
                         qualifiers={'locus_tag': ['test'], **qualifiers})
    record.features = [feature]
    return record, feature


class TranslationTests(unittest.TestCase):
    def test_forward_reverse_frame_and_alternative_genetic_code(self):
        cases = [
            ('ATGAAATAA', 1, {}, 'MK'),
            ('CTATTTCAT', -1, {}, 'MK'),
            ('AATGAAATAA', 1, {'codon_start': ['2']}, 'MK'),
            ('TTGCCCTAG', 1, {'transl_table': ['11']}, 'MP'),
            ('ATGTGATAA', 1, {'transl_table': ['4']}, 'MW'),
        ]
        for dna, strand, qualifiers, expected in cases:
            with self.subTest(dna=dna, strand=strand, qualifiers=qualifiers):
                record, feature = feature_record(dna, strand, **qualifiers)
                self.assertEqual(runtime.recover_translation(record, feature), (expected, ''))

    def test_absent_translation_table_uses_genbank_standard_code(self):
        record, feature = feature_record('ATGTGATAA')
        protein, reason = runtime.recover_translation(record, feature)
        self.assertIsNone(protein)
        self.assertIn('stop', reason.lower())

    def test_disrupted_ambiguous_or_exceptional_cds_is_not_translated(self):
        cases = [
            ('ATGAAATAA', {'pseudo': ['']}, 'pseudogene'),
            ('ATGAAATAA', {'pseudogene': ['unitary']}, 'pseudogene'),
            ('ATGAAATAA', {'transl_except': ['(pos:4..6,aa:Sec)']}, 'translation_exception'),
            ('ATGAAATAA', {'exception': ['RNA editing']}, 'translation_exception'),
            ('ATGAAATAA', {'note': ['frameshifted gene']}, 'noted_incomplete_or_disrupted'),
            ('ATGTAGAAATAA', {}, 'invalid_cds:'),
            ('ATGAAA', {}, 'invalid_cds:'),
            ('ATGAAATAAA', {}, 'invalid_cds:'),
            ('AAATAA', {}, 'invalid_cds:'),
            ('ATGNNNTAA', {}, 'missing_or_ambiguous_nucleotides'),
            ('ATGAAATAA', {'transl_table': ['9999']}, 'invalid_cds:'),
            ('ATGAAATAA', {'codon_start': ['4']}, 'invalid_cds:'),
        ]
        for dna, qualifiers, expected in cases:
            with self.subTest(dna=dna, qualifiers=qualifiers):
                record, feature = feature_record(dna, **qualifiers)
                protein, reason = runtime.recover_translation(record, feature)
                self.assertIsNone(protein)
                self.assertTrue(reason.startswith(expected), reason)

    def test_partial_compound_remote_and_out_of_bounds_calls_are_excluded(self):
        locations = [
            (SimpleLocation(BeforePosition(0), 9, strand=1), 'partial_or_uncertain_location'),
            (CompoundLocation([SimpleLocation(0, 3, strand=1), SimpleLocation(6, 9, strand=1)]), 'compound_or_remote_location'),
            (SimpleLocation(0, 9, strand=1, ref='elsewhere'), 'compound_or_remote_location'),
            (SimpleLocation(0, 12, strand=1), 'location_outside_sequence'),
            (SimpleLocation(0, 9), 'unknown_strand'),
            (None, 'missing_location'),
        ]
        for location, expected in locations:
            with self.subTest(location=location):
                record, feature = feature_record()
                feature.location = location
                self.assertEqual(runtime.recover_translation(record, feature), (None, expected))

    def test_record_without_nucleotide_sequence_cannot_be_recovered(self):
        record, feature = feature_record()
        record.seq = Seq(None, length=9)
        protein, reason = runtime.recover_translation(record, feature)
        self.assertIsNone(protein)
        self.assertTrue(reason.startswith('invalid_cds:'), reason)

    def test_all_missing_translations_recovered_and_original_upload_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); inputs = root / 'input'; inputs.mkdir()
            records = [feature_record(), feature_record('CTATTTCAT', -1),
                       feature_record('TTGCCCTAG', transl_table=['11'])]
            original = inputs / 'bakta.gb'
            SeqIO.write([record for record, _ in records], original, 'genbank')
            original_bytes = original.read_bytes()
            extension, paths = worker.prepare_bins(inputs, root / 'bins')
            self.assertEqual(runtime.prepare_genbank(root / 'bins', extension, root / 'out'), [])
            self.assertEqual(original.read_bytes(), original_bytes)
            with paths[0].open() as handle:
                prepared = list(SeqIO.parse(handle, 'genbank'))
            self.assertEqual([r.features[0].qualifiers['translation'][0] for r in prepared], ['MK', 'MK', 'MP'])
            with (root / 'out/genbank-translations.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[-1]['translation_table'], '11')
            self.assertEqual(rows[-1]['protein_length'], '2')
            self.assertTrue(all(row['source'] == 'annotated_CDS_translation' for row in rows))

    def test_supplied_translation_is_kept_even_when_dna_differs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); bins = root / 'bins'; bins.mkdir()
            record, _ = feature_record(translation=['MORIGINAL'])
            genome = bins / 'mpao1.gbk'; SeqIO.write(record, genome, 'genbank')
            original = genome.read_bytes()
            runtime.prepare_genbank(bins, 'gbk', root / 'out')
            self.assertEqual(genome.read_bytes(), original)
            self.assertEqual(SeqIO.read(genome, 'genbank').features[0].qualifiers['translation'], ['MORIGINAL'])

    def test_exclusion_report_is_written_when_no_proteins_can_be_recovered(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); bins = root / 'bins'; bins.mkdir()
            record, _ = feature_record(pseudogene=['unitary'])
            SeqIO.write(record, bins / 'sample.gbk', 'genbank')
            with self.assertRaisesRegex(RuntimeError, 'no usable CDS'):
                runtime.prepare_genbank(bins, 'gbk', root / 'out')
            with (root / 'out/genbank-exclusions.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]['reason'], 'pseudogene')


if __name__ == '__main__':
    unittest.main()
