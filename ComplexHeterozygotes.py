import re
from enum import Enum


class Ancestry(Enum):
    EA1KG = 1
    AAEVS = 2
    AF = 3
    EUEVS = 4
    GMAF = 5
    AFR1KGAF = 6

    @staticmethod
    def all():
        return [member for name, member in Ancestry.__members__.items()]


class Header:
    def __init__(self, header_line):
        self.titles = header_line.split("\t")
        self.patient_match = re.compile("[A-Z]{2}\d{5}[.][A-Z]{2}")
        self.patient_columns = [x for x in self.titles if self.patient_match.match(x) is not None]

    def __getitem__(self, item):
        return self.titles[item]

    def __len__(self):
        return len(self.titles)

    def __repr__(self):
        return "\t".join(self.titles) + "\n"


class Mutation:
    def __init__(self, header, line):
        columns = line.split("\t")
        assert len(columns) == len(header)
        self.data = {}
        for i in range(len(columns)):
            self.data[header[i]] = columns[i]
        self.header = header

    def var_recode(self):
        """
        This function recodes all the patient's alleles, to the number of variant alleles
        ex; if 'C' is reference
            "CC" -> '0'
            "AC" -> '1'
            "AA" -> '2'
        """
        var = self.data["var_allele"]
        for patient in self.header.patient_columns:
            alleles = self.data[patient]
            if alleles == "NA":
                self.data[patient] = "NA"
            else:
                if self.data["muttype"] == "SNP":
                    number_variant_chromosomes = alleles.count(var)
                    self.data[patient] = number_variant_chromosomes
                elif self.data["muttype"] == "INDEL":
                    number_variant_chromosomes = alleles.split(':').count(var)
                    self.data[patient] = number_variant_chromosomes
                else:
                    raise IOError("Don't know what to do with %s" % self.data["type"])

    def add_frequencies(self):
        """
        Adds total statistics
        """
        self.data['AC'] = 0
        self.data['AMax'] = 0
        for patient in self.header.patient_columns:
            try:
                self.data['AC'] += int(self.data[patient])
            except ValueError:
                pass
            self.data['AMax'] += 2
        self.data['AF'] = float(self.data['AC']) / float(self.data['AMax'])

    def rare_variant(self, threshold=0.02, populations_to_consider=Ancestry.all()):
        """
        :populations_to_consider: A list of objects of type Ancestry
        :return:Boolean, if this mutation is "rare"; <.02 for all populations
        """
        if Ancestry.EA1KG in populations_to_consider:
            try:
                allele_frequency = float(self.data["Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)"])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except ValueError:
                pass
        if Ancestry.AAEVS in populations_to_consider:
            try:
                allele_frequency = float(
                    self.data["Minor Allele Frequency in percent in the order of EA,AA,All"].split(",")[1])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.AF in populations_to_consider:
            pass
        if Ancestry.EUEVS in populations_to_consider:
            try:
                allele_frequency = float(
                    self.data["Minor Allele Frequency in percent in the order of EA,AA,All"].split(",")[0])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.GMAF in populations_to_consider:
            try:
                allele_frequency = float(
                    self.data["Estimated allele frequency in the range (0,1)"]) / 100
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.AFR1KGAF in populations_to_consider:
            try:
                allele_frequency = float(
                    self.data["Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)"]) / 100
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        return True

    def __repr__(self):
        return '\t'.join([str(self.data[title]) for title in self.header.titles]) + "\n"


class Data:
    def __init__(self, annotated_file_path):
        self.annotated_file_path = annotated_file_path
        self.mutations = []
        self.header = None
        self.load_mutations()
        self.var_recode()
        self.add_frequencies()
        self.remove_common_variants()

    def load_mutations(self):
        self.mutations = []
        with open(self.annotated_file_path) as f:
            self.header = Header(f.readline()[:-1])
            for line in f:
                self.mutations.append(Mutation(self.header, line[:-1]))

    def var_recode(self):
        for mutation in self.mutations:
            mutation.var_recode()

    def add_frequencies(self):
        self.header.titles.append("AC")
        self.header.titles.append("AMax")
        self.header.titles.append("AF")
        for mutation in self.mutations:
            mutation.add_frequencies()

    def remove_common_variants(self):
        self.mutations = [mutation for mutation in self.mutations if mutation.rare_variant()]

    def save_to_txt(self, output_file_path):
        with open(output_file_path, 'w') as f:
            f.write(self.header.__repr__())
            for mutation in self.mutations:
                f.write(mutation.__repr__())

if __name__ == '__main__':
    Data(r'AnnotatedTest.txt').save_to_txt('OutputTest.txt')
