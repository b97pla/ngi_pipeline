import locale
import StringIO
import unittest

from ngi_pipeline.engines.sarek.parsers import QualiMapParser, PicardMarkDuplicatesParser


class TestQualiMapParser(unittest.TestCase):

    def setUp(self):
        locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        self.examples = {
            "numeric_assignment": [
                "     std insert size = 1,478,053.3171",
                "     number of sequenced bases = 131,169,423 bp",
                "     number of reads = 967,817,737",
                "     GC percentage = 47.15%"],
            "assignment": [
                "     outfile = U004/genome_results.txt"],
            "cumulative_coverage": [
                "     There is a 0% of reference with a coverageData >= 49X",
                "     There is a 1.07% of reference with a coverageData >= 1X"],
            "contig_coverage": [
                "        chr1    248956422       58788   2.3613771248688655E-4   0.029014388002610237",
                "        chr1_KI270706v1_random  175055  0       0.0     0.0",
                "        21    77554325       780800934   10.067793562770355   1.07"],
            "invalid": [
                "",
                "                     ",
                ">>>>>>> Coverage per contig"
            ]
        }
        self.examples_expected = {
            "numeric_assignment": [
                {"std insert size": 1478053.3171},
                {"number of sequenced bases": 131169423},
                {"number of reads": 967817737},
                {"GC percentage": 47.15}],
            "assignment": [
                {"outfile": "U004/genome_results.txt"}],
            "cumulative_coverage": [
                {"49X": 0},
                {"1X": 1.07}],
            "contig_coverage": [
                {"chr1 coverage": {
                    "contig": "chr1",
                    "length": 248956422,
                    "mapped bases": 58788,
                    "mean coverage": 2.3613771248688655E-4,
                    "standard deviation": 0.029014388002610237}},
                {"chr1_KI270706v1_random coverage": {
                    "contig": "chr1_KI270706v1_random",
                    "length": 175055,
                    "mapped bases": 0,
                    "mean coverage": 0,
                    "standard deviation": 0}},
                {"21 coverage": {
                    "contig": "21",
                    "length": 77554325,
                    "mapped bases": 780800934,
                    "mean coverage": 10.067793562770355,
                    "standard deviation": 1.07}}]}

    def tearDown(self):
        locale.resetlocale()

    def helper(self, test_key, test_fn):
        for i in range(len(self.examples[test_key])):
            self.assertDictEqual(
                self.examples_expected[test_key][i],
                test_fn(self.examples[test_key][i]) or dict())
        for key in filter(lambda k: k != test_key, self.examples.keys()):
            self.assertTrue(
                all(map(
                    lambda l: test_fn(l) is None,
                    self.examples[key])))

    def test___parse_numeric_assignment(self):
        test_key = "numeric_assignment"
        test_fn = QualiMapParser._parse_numeric_assignment
        self.helper(test_key, test_fn)

    def test___parse_assignment(self):
        test_key = "assignment"
        test_fn = QualiMapParser._parse_assignment
        for i in range(len(self.examples[test_key])):
            self.assertDictEqual(
                self.examples_expected[test_key][i],
                test_fn(self.examples[test_key][i]) or dict())

    def test___parse_cumulative_coverage(self):
        test_key = "cumulative_coverage"
        test_fn = QualiMapParser._parse_cumulative_coverage
        self.helper(test_key, test_fn)

    def test___parse_contig_coverage(self):
        test_key = "contig_coverage"
        test_fn = QualiMapParser._parse_contig_coverage
        self.helper(test_key, test_fn)

    def test_get_autosomal_coverage(self):
        parser = QualiMapParser()
        QualiMapParser.AUTOSOMES = ["1", "21"]
        for contig in self.examples_expected["contig_coverage"]:
            parser.data.update(contig)
        expected_coverage = float(
            1. *
            sum([parser.data["{} coverage".format(c)]["mapped bases"] for c in ["chr1", "21"]]) /
            sum([parser.data["{} coverage".format(c)]["length"] for c in ["chr1", "21"]]))
        self.assertEqual(expected_coverage, parser.get_autosomal_coverage())

    def test_get_total_reads(self):
        expected_number_of_reads = 1234567890
        parser = QualiMapParser()
        parser.data["number of reads"] = expected_number_of_reads
        self.assertEqual(expected_number_of_reads, parser.get_total_reads())


class TestPicardMarkDuplicatesParser(unittest.TestCase):

    def setUp(self):
        self.example_output = """
        ## htsjdk.samtools.metrics.StringHeader
        # picard.sam.MarkDuplicates INPUT=[/scratch/ANALYSIS/EQ-5678/piper_ngi/05_processed
        ## htsjdk.samtools.metrics.StringHeader
        # Started on: Tue Sep 04 13:35:35 CEST 2018
        
        ## METRICS CLASS        picard.sam.DuplicationMetrics
        LIBRARY UNPAIRED_READS_EXAMINED READ_PAIRS_EXAMINED     UNMAPPED_READS  UNPAIRED_READ_DUPLICATES        READ_PAIR_DUPLICATES    READ_PAIR_OPTICAL_DUPLICATES    PERCENT_DUPLICATION     ESTIMATED_LIBRARY_SIZE
        SXYX_default-SeqLib.v22 15      1454662 19      0       62237   5723    0.042784        18088186
        SXYX_default-SeqLib.v23 12      484880  14      0       7033    1845    0.014504        22325476
        """
        self.example_handle = StringIO.StringIO(self.example_output)
        self.example_data = [
            {
                "LIBRARY": "SXYX_default-SeqLib.v22",
                "UNPAIRED_READS_EXAMINED": 15,
                "READ_PAIRS_EXAMINED": 1454662,
                "UNMAPPED_READS": 19,
                "UNPAIRED_READ_DUPLICATES": 0,
                "READ_PAIR_DUPLICATES": 62237,
                "READ_PAIR_OPTICAL_DUPLICATES": 5723,
                "PERCENT_DUPLICATION": 0.042784,
                "ESTIMATED_LIBRARY_SIZE": 18088186
            },
            {
                "LIBRARY": "SXYX_default-SeqLib.v23",
                "UNPAIRED_READS_EXAMINED": 12,
                "READ_PAIRS_EXAMINED": 484880,
                "UNMAPPED_READS": 14,
                "UNPAIRED_READ_DUPLICATES": 0,
                "READ_PAIR_DUPLICATES": 7033,
                "READ_PAIR_OPTICAL_DUPLICATES": 1845,
                "PERCENT_DUPLICATION": 0.014504,
                "ESTIMATED_LIBRARY_SIZE": 22325476
            }]

    def test__parse_metrics_handle(self):
        parser = PicardMarkDuplicatesParser()
        parser._parse_metrics_handle(self.example_handle)
        number_of_example_libs = len(self.example_data)
        self.assertEqual(number_of_example_libs, len(parser.data))
        for i in range(number_of_example_libs):
            self.assertDictEqual(self.example_data[i], parser.data[i])

    def test_get_percent_duplication(self):
        parser = PicardMarkDuplicatesParser()
        parser.data = self.example_data
        epsilon = 1e-4
        self.assertLess(
            abs(
                100*self.example_data[1]["PERCENT_DUPLICATION"] -
                parser.get_percent_duplication(library=self.example_data[1]["LIBRARY"])),
            epsilon)
        total_duplication = 3.571436857568654
        self.assertLess(
            abs(
                total_duplication -
                parser.get_percent_duplication()
            ),
            epsilon
        )

    def test__convert_to_unit(self):
        test_data = {
            int: [
                ["5", 5],
                ["-5", -5]],
            float: [
                ["3.5", 3.5],
                ["-3.5", -3.5]],
            str: [
                ["not-a-value", "not-a-value"],
                ["3,5", "3,5"]]}

        for expected_type in test_data.keys():
            for case_data in test_data[expected_type]:
                observed_value = PicardMarkDuplicatesParser._convert_to_unit(case_data[0])
                self.assertTrue(isinstance(observed_value, expected_type))
                self.assertEqual(case_data[1], observed_value)
