import re
from enum import Enum
import os


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
        self.patient_match = re.compile("[A-Z]{2}\d{5}[.][A-Z]{2}")  # Matches titles that look like 'HG00590.NA'
        self.patient_columns = [x for x in self.titles if self.patient_match.match(x) is not None]
        self.count_columns = ["Coding_gene"] + self.patient_columns
        self.var_columns = [x for x in self.titles if x not in self.count_columns] + ["AC", "AMax", "AF"]

    def __getitem__(self, item):
        return self.titles[item]

    def __len__(self):
        return len(self.titles)

    def __repr__(self):
        return "\t".join(self.titles) + os.linesep


class Mutation:
    def __init__(self, header, line, abbreviated_titles):
        columns = line.split("\t")
        assert len(columns) == len(header)
        self.data = {}
        for i in range(len(columns)):
            self[header[i]] = columns[i]
        self.header = header
        self.abbreviated_titles = abbreviated_titles

    def var_recode(self):
        """
        This function recodes all the patient's alleles, to the number of variant alleles
        ex; if 'C' is reference and 'A' is variant
            "CC" -> '0'
            "AC" -> '1'
            "AA" -> '2'
        If there are more variants than reference alleles, this will swap them
        ex; if 'C' is reference and 'A' is variant, but almost every patient has 'AA'
        This will make 'A' reference and 'C' variant, and count them that way.
        This is better than Allison's previous method, because it correctly counts other mutations
        ex; 'G' will still be labeled as a variant, regardless of whether 'C' or 'A' is reference
        """
        ref = self["ref_allele"]
        var = self["var_allele"]
        counts = {'ref': 0, 'var': 0}
        flip_data = dict(self.data)  # If the reference and variant alleles were mixed up
        flip_data["ref_allele"] = self["var_allele"]
        flip_data["var_allele"] = self["ref_allele"]
        for patient in self.header.patient_columns:
            alleles = self[patient]
            if alleles is not "NA":
                if self["muttype"] == "SNP":
                    number_variant_chromosomes = alleles.count(var)
                    number_reference_chromosomes = alleles.count(ref)
                elif self["muttype"] == "INDEL":
                    number_variant_chromosomes = alleles.split(':').count(var)
                    number_reference_chromosomes = alleles.split(':').count(ref)
                else:
                    raise IOError("Don't know what to do with %s" % self["type"])
                self.data[patient] = number_variant_chromosomes
                flip_data[patient] = number_reference_chromosomes
                counts['ref'] += number_reference_chromosomes
                counts['var'] += number_variant_chromosomes
        if counts['var'] > counts['ref']:
            self.data = flip_data

    def add_frequencies(self):
        """
        Adds total statistics: Number of variants, number of reference, frequency of variants
        """
        self['AC'] = 0
        self['AMax'] = 0
        for patient in self.header.patient_columns:
            try:
                self['AC'] += int(self[patient])
            except ValueError:
                pass
            self['AMax'] += 2
        self['AF'] = float(self['AC']) / float(self['AMax'])

    def rare_variant(self, threshold=0.02, populations_to_consider=Ancestry.all()):
        """
        :populations_to_consider: A list of objects of type Ancestry
        :return: Boolean, if this mutation is "rare"; <.02 for all populations
        """
        if Ancestry.EA1KG in populations_to_consider:
            try:
                allele_frequency = float(self["Allele frequency in the EUR"])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except ValueError:
                pass
        if Ancestry.AAEVS in populations_to_consider:
            try:
                allele_frequency = float(
                    self["Minor Allele Frequency in percent"].split(",")[1])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.AF in populations_to_consider:
            pass
        if Ancestry.EUEVS in populations_to_consider:
            try:
                allele_frequency = float(
                    self["Minor Allele Frequency in percent"].split(",")[0])
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.GMAF in populations_to_consider:
            try:
                allele_frequency = float(
                    self["Estimated allele frequency"]) / 100
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        if Ancestry.AFR1KGAF in populations_to_consider:
            try:
                allele_frequency = float(
                    self["Allele frequency in the AFR populations"]) / 100
                if threshold < allele_frequency < (1 - threshold):
                    return False
            except (ValueError, IndexError):
                pass
        return True

    def which_gene(self):
        effects = self["Predicted effects for this variant"].split(",")
        for effect in effects:
            summary = effect[:effect.index("(")]
            if summary in ['CODON_CHANGE_PLUS_CODON_DELETION', 'CODON_CHANGE_PLUS_CODON_INSERTION', 'CODON_INSERTION',
                           'CODON_DELETION', 'NON_SYNONYMOUS_CODING', 'STOP_GAINED', 'STOP_LOST', 'START_GAINED',
                           'START_LOST', 'FRAME_SHIFT', 'SPLICE_SITE_ACCEPTOR', 'SPLICE_SITE_DONOR',
                           'NON_SYNONYMOUS_START', 'EXON_DELETION']:
                self["Coding_gene"] = effect.split('|')[5]
                return

    def __repr__(self):
        return '\t'.join([str(self.data[title]) if title in self.data.keys() else "" for title in self.header.titles]) + os.linesep

    def __getitem__(self, item):
        """
        This takes advantage of the complicated abbreviated title lookup.
        If you're searching for item...
        -Check if item is in the dictionary
            -If so, return it
            -If not, check if item is in the abbreviated lookup dictionary (ie, 'Ancestry' for 'Ancestry Alleles')
                -If so, lookup the full key, then use that to lookup the value
                -If not, look to see if item is a new, previously unseen, abbreviation
                    -If so, put a new entry in the abbreviated lookup dictionary
                        (ie self.abbreviated_titles['Ancestry'] = 'Ancestry Alleles')
                    -If item matches too many or too few entries, raise an Exception
        """
        if item in self.data:
            return self.data[item]
        else:
            if item in self.abbreviated_titles:
                return self.data[self.abbreviated_titles[item]]
            else:
                possible_key_matches = []
                for key in self.data.keys():
                    if item in key:
                        possible_key_matches.append(key)
                if len(possible_key_matches) == 1:  # This is a substring of one and only one column
                    self.abbreviated_titles[item] = possible_key_matches[0]
                    return self.data[possible_key_matches[0]]
                elif len(possible_key_matches) == 0:
                    raise KeyError("No matches found.  Mispelling?")
                else:
                    raise KeyError("Matches more than one entry")

    def __setitem__(self, key, value):
        self.data[key] = value


class CHCount(Mutation):
    def __init__(self, ch_count_header, *args):
        self.header = ch_count_header
        if len(args) == 1:
            if isinstance(args[0], dict):
                self.data = args[0]
                return
        elif len(args) == 2:
            if isinstance(args[0], Header) and isinstance(args[1], Mutation):
                data_header, mutation = args
                self.data = {key: mutation.data[key] for key in data_header.count_columns}
                return
        print args
        raise TypeError


class CHCountHeader(Header):
    def __init__(self, header):
        self.titles = header.count_columns


class CHVariant(Mutation):
    def __init__(self, ch_variant_header, data_header, mutation):
        self.data = {key: mutation.data[key] for key in data_header.var_columns}
        self.header = ch_variant_header


class CHVariantHeader(Header):
    def __init__(self, header):
        self.titles = header.var_columns

class Data:
    def __init__(self, annotated_file_path):
        self.abbreviated_titles = {}
        """
        Abbreviated titles is a dictionary that is shallow copied into all the instances of Mutation
        It's role is to provide an easy lookup for abbreviated titles that is shared between mutations
        The first you look up mutation['Ancestral'] instead of mutation["Ancestral Allele"], the program will
            have to look through all the keys to find if any match
        Afterwards, however, an entry will be added to abbreviated_titles;
            abbreviated_titles['Ancestral'] = 'Ancestral Allele'
        Henceforth, the program will use this lookup if possible, instead of having to iterate through all the keys
        """
        self.annotated_file_path = annotated_file_path
        self.mutations = []
        self.header = None

        self.ch_counts = []
        self.ch_count_header = None

        self.ch_variants = []
        self.ch_variant_header = None

        self.load_mutations()
        self.var_recode()
        self.add_frequencies()
        self.remove_common_variants()
        self.remove_repeats()
        self.which_gene()
        self.convert_to_ch()
        self.concatenate_counts()

    def load_mutations(self):
        """
        This function loads the original input file of the program
        The header (titles) will be put into self.header
        Each line, detailing one mutation, will be put into a Mutation class, and placed in self.mutations

        Mutation class is a wrapper around a dictionary, and Header is a wrapper around a list.
        Header has the list of titles, in order.
        Mutation's dictionary is intrinsically unordered.
        To write the first column;
            The title is the first entry in Header's list
            The column values can be looked up in Mutation's dictionary
        """
        self.mutations = []
        with open(self.annotated_file_path) as f:
            self.header = Header(f.readline()[:-1])
            for line in f:
                self.mutations.append(Mutation(self.header, line[:-1], self.abbreviated_titles))

    def var_recode(self):
        """
        Shadows varRecode1.py from Allison's scripts
        This replaces the actual alleles ('AA') with the number of variant alleles ('2')
        """
        for mutation in self.mutations:
            mutation.var_recode()

    def add_frequencies(self):
        """
        This function computes some simple statistics on each mutation.
        Computes the number of mutated alleles, the number of wild-type alleles, and the frequency of mutations
        """
        self.header.titles.append("AC")
        self.header.titles.append("AMax")
        self.header.titles.append("AF")
        for mutation in self.mutations:
            mutation.add_frequencies()

    def remove_common_variants(self):
        """
        Removes all "common" variants
        that have an incidence of 2% or greater in the populations
        Currently, this function excludes variants that are 2% or greater in any of the populations.
        However, mutation.rare_variant takes an optional argument which can specify which populations should be checked
        """
        self.mutations = [mutation for mutation in self.mutations if mutation.rare_variant()]

    def remove_repeats(self):
        """
        Code to remove duplicates from
        http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
        Removes mutations with the same "dbID"
        Currently, this just deletes the second instance.
        TODO: Ask Allison if the mutations should be concatenated somehow
        """
        seen = set()
        unique_mutations = []
        seen_add = seen.add
        for mutation in self.mutations:
            if mutation['dbID'] not in seen:
                seen_add(mutation['dbID'])  # Python is stupid, and will check the object each time otherwise
                unique_mutations.append(mutation)
        print "%s mutations, %s unique" % (len(self.mutations), len(unique_mutations))
        self.mutations = unique_mutations

    def which_gene(self):
        """
        Adds the "Coding_gene" column.
        """
        self.header.titles.insert(4, "Coding_gene")
        for mutation in self.mutations:
            mutation.which_gene()
        self.mutations.sort(key=lambda x: x["Coding_gene"])

    def convert_to_ch(self):
        self.ch_counts = []
        self.ch_variants = []
        self.ch_count_header = CHCountHeader(self.header)
        self.ch_variant_header = CHVariantHeader(self.header)
        for mutation in self.mutations:
            self.ch_counts.append(CHCount(self.ch_count_header, self.header, mutation))
            self.ch_variants.append(CHVariant(self.ch_variant_header, self.header, mutation))

    def concatenate_counts(self):
        self.ch_counts.sort(key=lambda x: x['Coding_gene'])
        begin_sequence = 0
        concatenated_counts = []
        while begin_sequence < len(self.ch_counts):
            end_sequence = begin_sequence + 1
            while self.ch_counts[end_sequence]["Coding_gene"] == self.ch_counts[begin_sequence]["Coding_gene"]:
                end_sequence += 1
                if end_sequence > len(self.ch_counts) - 1:
                    break
            concatenated_data = {"GENE": self.ch_counts[begin_sequence]["Coding_gene"], "CHCOUNT": 0}
            for patient in self.header.patient_columns:
                concatenated_data[patient] = 0
                for index in range(begin_sequence, end_sequence):
                    if int(self.ch_counts[index][patient]) > concatenated_data[patient]:
                        concatenated_data[patient] = int(self.ch_counts[index][patient])
                    if self.ch_counts[index][patient] != 0:
                        concatenated_data["CHCOUNT"] += 1
            concatenated_counts.append(CHCount(self.ch_count_header, concatenated_data))
            begin_sequence = end_sequence
        self.ch_counts = concatenated_counts
        self.ch_count_header.titles[0] = "GENE"
        self.ch_count_header.titles.insert(1, "CHCOUNT")

    def save_to_txt(self, output_file_path):
        """
        Write the mutations and the headers to output_file_path
        """
        with open(output_file_path, 'w') as f:
            f.write(self.header.__repr__())
            for mutation in self.mutations:
                f.write(mutation.__repr__())

    def save_ch_files(self, var_file_path, count_file_path):
        with open(var_file_path, 'w') as f:
            f.write(self.ch_variant_header.__repr__())
            for ch_variant in self.ch_variants:
                f.write(ch_variant.__repr__())
        with open(count_file_path, 'w') as f:
            f.write(self.ch_count_header.__repr__())
            for ch_count in self.ch_counts:
                f.write(ch_count.__repr__())

if __name__ == '__main__':
    Data(r'AnnotatedTest.txt').save_ch_files('Var_test.txt', 'Count_test.txt')
