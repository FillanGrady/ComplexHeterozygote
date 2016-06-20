import re
from enum import Enum
import os
import time
from multiprocessing.dummy import Pool
import functools


class Ancestry(Enum):
    Overall = 1
    African = 2
    American = 3
    EastAsian = 4
    European = 5
    SouthAsian = 6

    @staticmethod
    def all():
        return [member for name, member in Ancestry.__members__.items()]


class Genotypes(Enum):
    """
    Exact comparisons only;
    OneMutation is for patients with one and only one mutation
    """
    CompoundHeterozygotes = 1
    Homozygous = 2
    TwoMutations = 3
    OneMutation = 4
    WildType = 5


class Gene(Enum):
    """
    What genes are affected by a mutation
    """
    INTRON = 1
    EXON = 2
    DOWNSTREAM = 3
    UPSTREAM = 4
    SYNONYMOUS_CODING = 5
    NON_SYNONYMOUS_CODING = 6
    UTR_3_PRIME = 7
    UTR_5_PRIME = 8
    START_GAINED = 9
    START_LOST = 10
    STOP_GAINED = 11
    CODON_DELETION = 12
    CODON_CHANGE_PLUS_CODON_DELETION = 13
    CODON_INSERTION = 14
    FRAME_SHIFT = 15
    SPLICE_SITE_DONOR = 16
    SPLICE_SITE_ACCEPTOR = 17

    @staticmethod
    def all():
        return [member for name, member in Gene.__members__.items()]

    @staticmethod
    def default():
        return [Gene.EXON, Gene.NON_SYNONYMOUS_CODING]

    def __str__(self):
        return self.name


class Header:
    def __init__(self, header_line):
        self.titles = header_line.split("\t")
        self.patient_match = re.compile("[A-Z]{2}\d{5}")  # Matches titles that look like 'HG00590'
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
    def __init__(self, *args):
        if len(args) == 3:
            header, line, abbreviated_titles = args
            if isinstance(header, Header) and isinstance(line, str) and isinstance(abbreviated_titles, dict):
                columns = line.split("\t")
                assert len(columns) == len(header)
                self.data = {}
                for title, column in zip(header.titles, columns):
                    self[title] = column
                self.header = header
                self.abbreviated_titles = abbreviated_titles
            else:
                raise ValueError("Unrecognized arguments")
        elif len(args) == 1:
            mutation = args[0]
            if isinstance(mutation, Mutation):
                self.header = mutation.header
                self.abbreviated_titles = mutation.abbreviated_titles
                self.data = dict(mutation.data)
            else:
                raise ValueError("Unrecognized arguments")
        else:
            raise ValueError("Wrong number of arguments")

    def add_frequencies(self):
        """
        Adds total statistics: Number of variants, number of reference, frequency of variants
        """
        self['AC'] = 0
        for patient in self.header.patient_columns:
            if self[patient][0] == '1':
                self['AC'] += 1
            if self[patient][2] == '1':
                self['AC'] += 1
        self['AMax'] = 2 * len(self.header.patient_columns)
        self['AF'] = float(self['AC']) / float(self['AMax'])
        """
        if self['AF'] > .5:  # Var and Ref need to be switched
            self['REF'], self['ALT'] = self['ALT'], self['REF']
            self['AC'] = self['AMax'] - self['AC']
            self['AF'] = 1 - self['AF']
            for patient in self.header.patient_columns:
                if self[patient] == '0|0':
                    self[patient] = '1|1'
                elif self[patient] == '1|0':
                    self[patient] = '0|1'
                elif self[patient] == '0|1':
                    self[patient] = '1|0'
                elif self[patient] == '1|1':
                    self[patient] = '0|0'
                else:
                    raise ValueError(self[patient])
        """

    def split_info(self):
        info = self["INFO"].split(";")
        del self.data["INFO"]
        for datum in info:
            datum = datum.strip('|').split("=")
            if len(datum) == 2:
                self[datum[0]] = datum[1]
                if datum[0] not in self.header.titles:
                    self.header.titles.insert(7, datum[0])

    def rare_variant(self, threshold=0.02, populations_to_consider=Ancestry.all()):
        """
        :populations_to_consider: A list of objects of type Ancestry
        :return: Boolean, if this mutation is "rare"; <.02 for all populations
        """
        try:
            if Ancestry.Overall in populations_to_consider:
                if threshold < float(self['AF']) < (1 - threshold):
                    return False
            if Ancestry.African in populations_to_consider:
                if threshold < float(self['AFR_AF']) < (1 - threshold):
                    return False
            if Ancestry.American in populations_to_consider:
                if threshold < float(self['AMR_AF']) < (1 - threshold):
                    return False
            if Ancestry.EastAsian in populations_to_consider:
                if threshold < float(self['EAS_AF']) < (1 - threshold):
                    return False
            if Ancestry.European in populations_to_consider:
                if threshold < float(self['EUR_AF']) < (1 - threshold):
                    return False
            if Ancestry.SouthAsian in populations_to_consider:
                if threshold < float(self['SAS_AF']) < (1 - threshold):
                    return False
        except ValueError:
            return False
        return True

    def get_coding_genes(self, affected_genes=Gene.default()):
        effects = self["EFF"].split(",")
        coding_genes = set()
        for effect in effects:
            summary = effect[:effect.index("(")]
            if summary in map(str, affected_genes):
                coding_genes.add(effect.split('|')[5])
        return coding_genes

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
        print(args)
        raise TypeError


class CHCountHeader(Header):
    def __init__(self, header):
        self.titles = header.count_columns
        self.patient_columns = header.patient_columns
        self.titles[0] = "Coding_gene"
        self.titles.insert(1, "CHCount")


class CHVariant(Mutation):
    def __init__(self, ch_variant_header, data_header, mutation):
        self.data = {key: mutation.data[key] for key in data_header.var_columns}
        self.header = ch_variant_header


class CHVariantHeader(Header):
    def __init__(self, header):
        self.titles = header.var_columns


class Data:
    def __init__(self, annotated_file_path):
        """
        Abbreviated titles is a dictionary that is shallow copied into all the instances of Mutation
        It's role is to provide an easy lookup for abbreviated titles that is shared between mutations
        The first you look up mutation['Ancestral'] instead of mutation["Ancestral Allele"], the program will
            have to look through all the keys to find if any match
        Afterwards, however, an entry will be added to abbreviated_titles;
            abbreviated_titles['Ancestral'] = 'Ancestral Allele'
        Henceforth, the program will use this lookup if possible, instead of having to iterate through all the keys
        """

        self.abbreviated_titles = {}
        self.annotated_file_path = annotated_file_path
        self.mutations = []
        self.header = None

        self.ch_counts = []
        self.ch_count_header = None

        start = time.time()
        self.load_mutations()
        print("Load mutations: %.2f" % (time.time() - start))
        self.add_frequencies()
        print("Add frequencies: %.2f" % (time.time() - start))
        self.split_info()
        print("Split info: %.2f" % (time.time() - start))
        self.remove_common_variants()
        print("Remove common variants: %.2f" % (time.time() - start))
        self.remove_repeats()
        print("Remove repeats: %.2f" % (time.time() - start))
        self.get_coding_genes()
        print("Get coding genes: %.2f" % (time.time() - start))
        gene_groups = [self.mutations[x[0]:x[1]] for x in self.group_genes()]
        print("Group genes: %.2f" % (time.time() - start))
        self.test_genotypes = [Genotypes.CompoundHeterozygotes]
        pool = Pool(4)
        self.ch_counts = pool.map(
            functools.partial(
                Data.concatenate_groups, ch_count_header=self.ch_count_header, test_genotypes=self.test_genotypes)
            , gene_groups)  # We can't lazily evaluate, because remove_var_patients will delete the records
        print("Concatenate groups: %.2f" % (time.time() - start))
        self.remove_var_patients()

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

    def split_info(self):
        """
        Column "INFO" looks like "AF=0.2;EUR_AF=0.3;EAS_AF=0.1..."
        split this into
        Column "AF" = 0.2
        Column "EUR_AF" = 0.3
        Column "EAS_AF" = 0.1
        ...
        """
        for mutation in self.mutations:
            mutation.split_info()

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
        """
        seen = set()
        unique_mutations = []
        seen_add = seen.add
        for mutation in self.mutations:
            if mutation['ID'] not in seen:
                seen_add(mutation['ID'])  # Python is stupid, and will check the object each time otherwise
                unique_mutations.append(mutation)
        self.mutations = unique_mutations

    def get_coding_genes(self):
        """
        Adds the "Coding_gene" column.
        """
        self.header.titles.insert(4, "Coding_gene")
        new_mutations = []
        for mutation in self.mutations:  # Iterate backwards through the list to avoid conflicts
            coding_genes = mutation.get_coding_genes()  # Sets remove duplicates from list
            for coding_gene in coding_genes:
                new_mutation = Mutation(mutation)  # TODO: To save memory, the dictionary doesn't need to be copied once
                new_mutation["Coding_gene"] = coding_gene
                new_mutations.append(new_mutation)
        self.mutations = sorted(new_mutations, key=lambda x: x['Coding_gene'])

    def group_genes(self):
        """
        Concatenates mutations that affect the same gene
        For example, if we have two mutations, that look like this
        ID      Coding_gene     Patient1    Patient2
        var1    gene1           1|0         1|1
        var2    gene2           0|1         0|0
        Our output would like this if Genotypes=[Genotypes.ComplexHeterozygous]
        Coding_gene     Patient1    Patient2
        gene1           1           0
        and this if Genotypes=[Genotypes.Homozygous]
        Coding_gene     Patient1    Patient2
        gene1           0           1

        This function simply groups mutations that affect the same gene.
        It'll return lists of tuples, where the first value is the start_index, and the second value is the end_index
        """
        self.ch_counts = []
        self.ch_count_header = CHCountHeader(self.header)
        self.mutations.sort(key=lambda x: x['Coding_gene'])
        gene_groups = []
        begin_sequence = 0
        end_sequence = 0
        while begin_sequence < len(self.mutations):
            while self.mutations[end_sequence]["Coding_gene"] == self.mutations[begin_sequence]["Coding_gene"]:
                end_sequence += 1
                if end_sequence > len(self.mutations) - 1:
                    break
            gene_groups.append((begin_sequence, end_sequence))
            begin_sequence = end_sequence
        return gene_groups

    @staticmethod
    def concatenate_groups(mutations, ch_count_header, test_genotypes):
        """
        Needs to be static for multiprocessing
        """
        concatenated_data = {'Coding_gene': mutations[0]["Coding_gene"], 'CHCount': 0}
        for patient in ch_count_header.patient_columns:
            maternal_mutations = set()
            paternal_mutations = set()
            for mutation_index, mutation in enumerate(mutations):
                alleles = mutation.data[patient].split('|')
                if alleles[0] != '0':
                    maternal_mutations.add(mutation_index)
                if alleles[1] != '0':
                    paternal_mutations.add(mutation_index)
            satisfies_any = False
            for test_genotype in test_genotypes:
                if test_genotype == Genotypes.CompoundHeterozygotes:
                    satisfies_any = satisfies_any or \
                                    ((len(maternal_mutations) > 0 and len(paternal_mutations) > 0) and
                                                          len(set.union(maternal_mutations, paternal_mutations)) >= 2)
                if test_genotype == Genotypes.Homozygous:
                    satisfies_any = satisfies_any or \
                                    len(set.intersection(maternal_mutations, paternal_mutations)) > 0
                if test_genotype == Genotypes.TwoMutations:
                    satisfies_any = satisfies_any or \
                                    (len(maternal_mutations) > 0 and len(paternal_mutations) > 0)
                if test_genotype == Genotypes.OneMutation:
                    satisfies_any = satisfies_any or \
                                    (len(maternal_mutations) > 0 != len(paternal_mutations) > 0)
                if test_genotype == Genotypes.WildType:
                    satisfies_any = satisfies_any or \
                                    (len(maternal_mutations) == 0 and len(paternal_mutations) == 0)
            concatenated_data[patient] = int(satisfies_any)
            concatenated_data["CHCount"] += concatenated_data[patient]
        return CHCount(ch_count_header, concatenated_data)

    def remove_var_patients(self):
        self.header.titles.insert(1, "Patients_affected")
        for mutation in self.mutations:
            patients_affected_by_mutation = []
            for patient in self.header.patient_columns:
                if mutation.data[patient] != "0|0":
                    patients_affected_by_mutation.append(patient)
                del mutation.data[patient]
            mutation["Patients_affected"] = ",".join(patients_affected_by_mutation)
        self.header.titles = [x for x in self.header.titles if x not in self.header.patient_columns]
        self.header.patient_columns = ["Patients_affected"]

    def save_to_txt(self, output_file_path):
        """
        Write the mutations and the headers to output_file_path
        """
        with open(output_file_path, 'w') as f:
            f.write(repr(self.header))
            for mutation in self.mutations:
                f.write(repr(mutation))

    def save_ch_files(self, count_file_path):
        with open(count_file_path, 'w') as f:
            f.write(repr(self.ch_count_header))
            for ch_count in self.ch_counts:
                f.write(repr(ch_count))

if __name__ == '__main__':
    d = Data(r'chr22.ANNOTATED.db.eff_coding.vcf')
    #d = Data(r'AnnotatedTest.vcf')
    start = time.time()
    d.save_to_txt("OutputTest.vcf")
    print("Save txt file: %.2f" % (time.time() - start))
    d.save_ch_files("CountTest.vcf")
    print("Save ch file: %.2f" % (time.time() - start))
