import locale
import unittest

from ngi_pipeline.engines.sarek.parsers import QualiMapParser


class TestQualiMapParser(unittest.TestCase):

    def setUp(self):
        locale.setlocale(locale.LC_ALL, 'en_US')
        self.examples = {
            "numeric_assignment": [
                "     std insert size = 1,478,053.3171",
                "     number of sequenced bases = 131,169,423 bp",
                "     GC percentage = 47.15%"],
            "assignment": [
                "     outfile = U004/genome_results.txt"],
            "cumulative_coverage": [
                "     There is a 0% of reference with a coverageData >= 49X",
                "     There is a 1.07% of reference with a coverageData >= 1X"],
            "contig_coverage": [
                "        chr1    248956422       58788   2.3613771248688655E-4   0.029014388002610237",
                "        chr1_KI270706v1_random  175055  0       0.0     0.0"],
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
                    "standard deviation": 0}}]}

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
