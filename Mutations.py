from Enumerations import Ancestry, GeneEffects, Genotypes
import os
import re


class Header:
    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], str):
                """
                This option reads the column headers in from an input file
                input_titles: read in from vcf file
                patient_columns: columns that match patient ids
                output_titles: to be output to a vcf file
                count_columns: to be output to the count file
                """
                header_line = args[0]
                self.input_titles = header_line.split("\t")
                info_columns = ["AA", "AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", "EAS_AF", "EUR_AF", "HRun", "NS",
                                "SAS_AF", "VT", "EFF"]
                self.patient_match = re.compile("[A-Z]{2}\d{5}")  # Matches titles that look like 'HG00590'
                self.patient_columns = [x for x in self.input_titles if self.patient_match.match(x) is not None]
                self.output_titles = self.input_titles[0:1] + ["PATIENTS"] + self.input_titles[1:4] + ["CODING_GENE"]\
                    + self.input_titles[4:7] + info_columns + \
                    [x for x in self.input_titles[7:] if x not in self.patient_columns] + \
                    ["AC", "AMax", "AF"]
                self.count_columns = ["CODING_GENE", "CHCOUNT"] + self.patient_columns
                self.var_columns = [x for x in self.output_titles if x not in self.count_columns]
            elif isinstance(args[0], Header):
                """
                This option creates a copy of an instance of the Header class
                """
                header = args[0]
                self.input_titles = header.input_titles
                self.output_titles = header.output_titles
                self.patient_columns = header.patient_columns
                self.count_columns = header.count_columns
                self.var_columns = header.var_columns
            else:
                raise ValueError("Unrecognized arguments")
        else:
            raise ValueError("Wrong number of arguments")

    def __getitem__(self, item):
        return self.output_titles[item]

    def __len__(self):
        return len(self.output_titles)

    def __repr__(self):
        return "\t".join(self.output_titles) + os.linesep


class Mutation:
    def __init__(self, *args):
        if len(args) == 3:
            """
            This option loads a line of the vcf file as a mutation
            The header (titles) have been put into self.header
            One line, detailing one mutation, will be put into this object, and placed in Data.mutation
            """
            header, line, abbreviated_titles = args
            if isinstance(header, Header) and isinstance(line, str) and isinstance(abbreviated_titles, dict):
                columns = line.split("\t")
                self.data = {}
                for title, column in zip(header.input_titles, columns):
                    self[title] = column
                self.header = header
                self.abbreviated_titles = abbreviated_titles
            else:
                raise ValueError("Unrecognized arguments")
        elif len(args) == 1:
            """
            This option create a copy of an instance of the Mutation class
            I couldn't figure out how to make __copy__ work
            """
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
        This function computes some simple statistics on each mutation.
        Computes the number of mutated alleles, the number of wild-type alleles, and the frequency of mutations
        """
        self['AC'] = 0
        for patient in self.header.patient_columns:
            if self[patient][0] == '1':
                self['AC'] += 1
            if self[patient][2] == '1':
                self['AC'] += 1
        self['AMax'] = 2 * len(self.header.patient_columns)
        self['AF'] = float(self['AC']) / float(self['AMax'])
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

    def split_info(self):
        """
        Column "INFO" looks like "AF=0.2;EUR_AF=0.3;EAS_AF=0.1..."
        split this into
        Column "AF" = 0.2
        Column "EUR_AF" = 0.3
        Column "EAS_AF" = 0.1
        ...
        """
        info = self["INFO"].split(";")
        del self.data["INFO"]
        for datum in info:
            datum = datum.strip('|').split("=")
            if len(datum) == 2:
                self[datum[0]] = datum[1]

    def rare_variant(self, threshold=0.02, populations_to_consider=Ancestry.all()):
        """
        This function returns whether or not this mutation is rare;
        if it has an incidence of 2% or greater in the populations
        However, mutation.rare_variant takes an optional argument which can specify which populations should be checked
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

    def get_coding_genes(self, affected_genes=GeneEffects.default()):
        effects = self["EFF"].split(",")
        coding_genes = set()
        for effect in effects:
            summary = effect[:effect.index("(")]
            if summary in map(str, affected_genes):
                coding_genes.add(effect.split('|')[5])
        return coding_genes

    def __repr__(self):
        return '\t'.join([str(self.data[title]) if title in self.data.keys() else ""
                          for title in self.header.output_titles]) + os.linesep

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
                    print(item)
                    raise KeyError("No matches found.  Misspelling?")
                else:
                    raise KeyError("Matches more than one entry")

    def __setitem__(self, key, value):
        self.data[key] = value


class PatientGenotype:
    """
    This class takes care of one coding_gene for one patient, for one list of genotypes
    For example, AP00025.1 for HG01302 and [ComplexHeterozygotes]
    Everytime you see a new mutation that will affect AP00025.1, call parse_new_alleles
    self.value holds whether or not this coding_gene for this patient is one of the types in valid_genotypes
    """
    def __init__(self, valid_genotypes):
        self.value = False  # Represents whether or not this patient has this genotype
        self.maternal_seen = False
        self.paternal_seen = False
        self.valid_genotypes = valid_genotypes

    def parse_new_alleles(self, alleles):
        """
        Update self.value with new information
        """
        if not self.value:
            maternal_allele = alleles[0] != "0"
            paternal_allele = alleles[1] != "0"  # False is wildtype, True is mutated
            if Genotypes.OneMutation in self.valid_genotypes:
                if maternal_allele or paternal_allele:
                    self.value = True
                    return
            if Genotypes.TwoMutations in self.valid_genotypes:
                if maternal_allele and (paternal_allele or self.paternal_seen):
                    self.value = True
                    return
                if paternal_allele and (maternal_allele or self.maternal_seen):
                    self.value = True
                    return
            if Genotypes.CompoundHeterozygotes in self.valid_genotypes:
                if maternal_allele and self.paternal_seen:
                    self.value = True
                    return
                if paternal_allele and self.maternal_seen:
                    self.value = True
                    return
            if Genotypes.Homozygous:
                if maternal_allele and paternal_allele:
                    self.value = True
                    return
            self.maternal_seen = self.maternal_seen or maternal_allele
            self.paternal_seen = self.paternal_seen or paternal_allele

    def __repr__(self):
        return str(int(self.value))


class CodingGeneMutation(Mutation):
    """
    This class holds genotype changes for every patient for one coding gene
    This class will have one PatientGenotype object for every patient,
    These PatientGenotype objects will be updated with every new mutation affecting the given coding_gene
    """
    def __init__(self, ch_count_header, coding_gene, valid_genotypes):
        self.header = ch_count_header
        self.valid_genotypes = valid_genotypes
        self.data = {'CODING_GENE': coding_gene, 'CHCOUNT': 0}
        for patient in self.header.patient_columns:
            self.data[patient] = PatientGenotype(self.valid_genotypes)

    def parse_new_mutation(self, mutation):
        for patient in self.header.patient_columns:
            self.data[patient].parse_new_alleles(mutation[patient])

    def __repr__(self):
        return '\t'.join(
            [str(self.data[title]) if title in self.data.keys() else "" for title in self.header.output_titles])\
               + os.linesep


class CodingGeneMutationHeader(Header):
    """
    Creates a header for the Counts
    """
    def __init__(self, header, my_patients):
        self.output_titles = header.count_columns
        self.patient_columns = header.patient_columns
        self.my_patients = my_patients

    def __repr__(self):
        patient_header = '\t'.join(self.output_titles)
        gender_header = '\t'.join([self.my_patients[patient].gender if patient in self.my_patients else ""
                                   for patient in self.output_titles])
        population_header = '\t'.join([self.my_patients[patient].population if patient in self.my_patients else ""
                                       for patient in self.output_titles])
        super_pop_header = '\t'.join([self.my_patients[patient].super_population if patient in self.my_patients else ""
                                      for patient in self.output_titles])
        return patient_header + os.linesep + gender_header + os.linesep + population_header + os.linesep + \
            super_pop_header
