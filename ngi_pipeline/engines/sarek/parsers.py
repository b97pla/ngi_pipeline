
import json
import locale
import re
import string


class MultiQCParser:

    def __init__(self):
        self.data = {}

    def data_source(self, tool, section):
        sources = self.data["report_data_sources"][tool][section]
        return list(sources.keys())

    def parse_json_data(self, json_data_file):
        with open(json_data_file) as fh:
            self.data = json.load(fh)

    def multiqc_raw_data(self):
        return self.data["report_saved_raw_data"]

    def multiqc_general_stats(self):
        return self.multiqc_raw_data()["multiqc_general_stats"]

    def multiqc_samtools_stats(self):
        return self.multiqc_raw_data()["multiqc_samtools_stats"]

    def autosomal_coverage(self):
        pass

    def total_sequenced_reads(self):
        source = list(filter(lambda s: s.endswith(".real"), self.data_source("Samtools", "stats"))).pop()
        source_data = self.multiqc_samtools_stats()[source]
        return source_data["raw_total_sequences"]


class QualiMapParser:

    def __init__(self):
        self.data = {}

    def parse_genome_results(self, genome_results_file):
        locale.setlocale(locale.LC_ALL, '')
        with open(genome_results_file) as fh:
            for line in fh:
                self.data.update(self._parse_entry(line))

    def get_autosomal_coverage(self):
        pass

    @staticmethod
    def _parse_numeric_assignment(line):
        # identify key = value assignments with numeric values
        try:
            key, value = re.search(r'^\s+([^=]+) = ([0-9,\.]+)', line).groups()
            value = locale.atof(value)
            key = key.strip()
            return {key: value}
        except ValueError:
            # it was an assignment but the conversion to numeric value failed
            pass
        except AttributeError:
            # not an assignment or not a numeric assignment
            pass

    @staticmethod
    def _parse_assignment(line):
        # identify key = value assignments
        try:
            key, value = re.search(r'^\s+([^=]+)= (.+)$', line).groups()
            key, value = map(string.strip, [key, value])
            return {key: value}
        except AttributeError:
            # not an assignment
            pass

    @staticmethod
    def _parse_cumulative_coverage(line):
        # identify cumulative coverage calculations
        try:
            value, key = re.search(r'([0-9\.]+)%.*>= ([0-9]+X)', line).groups()
            value = locale.atof(value)
            return {key: value}
        except AttributeError:
            pass

    @staticmethod
    def _parse_contig_coverage(line):
        # identify contig coverage
        try:
            values = re.search(r'^\s+(\S+)\s+([0-9]+)\s+([0-9]+)\s+(\S+)\s+(\S+)\s*$', line).groups()
            key = "{} coverage".format(values[0])
            value = dict()
            value["contig"] = values[0]
            value["length"] = locale.atoi(values[1])
            value["mapped bases"] = locale.atoi(values[2])
            value["mean coverage"] = locale.atof(values[3])
            value["standard deviation"] = locale.atof(values[4])
            return {key: value}
        except ValueError:
            # problems with the conversion to numeric values
            pass
        except AttributeError:
            # not a contig coverage row
            pass

    def _parse_entry(self, line):
        return self._parse_numeric_assignment(line) or \
               self._parse_assignment(line) or \
               self._parse_cumulative_coverage(line) or \
               self._parse_contig_coverage(line) or \
               dict()
