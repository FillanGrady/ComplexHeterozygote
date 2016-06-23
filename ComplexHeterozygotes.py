import sys
from Patients import Patients
from Enumerations import Genotypes
from Mutations import Mutation, Header, CodingGeneMutation, CodingGeneMutationHeader
import zipfile


class Data:
    """
    This class is a wrapper for a whole vcf file
    It reads in the input file line by line;
        It keeps track of the Header, and one Mutation at a time
    It also keeps a dictionary of coding_genes that have been affected by any of the mutations
        and constantly updates this dictionary as it sees new mutations
    It acts as a pipe for mutations, collecting their effects on coding_genes
    """
    def __init__(self, annotated_file_path, output_file_path, count_file_path, subject_info_file_path,
                 valid_genotypes=None):
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
        self.output_file_path = output_file_path
        self.subject_info_file_path = subject_info_file_path
        if valid_genotypes is None:
            self.valid_genotypes = [Genotypes.CompoundHeterozygotes]
        else:
            self.valid_genotypes = valid_genotypes
        self.mutation = None
        self.mutations = []
        self.header = None
        """
        Mutation class is a wrapper around a dictionary, and Header is a wrapper around a list.
        Header has the list of titles, in order. Mutation's dictionary is intrinsically unordered.
        To write the first column; The title is the first entry in Header's list
        The column values can be looked up in Mutation's dictionary
        """

        self.ch_counts = {}
        self.ch_count_header = None
        self.seen_ids = set()
        with open(self.output_file_path, 'w') as self.output_file:
            if zipfile.is_zipfile(self.annotated_file_path):
                with zipfile.ZipFile(self.annotated_file_path, mode='r') as zf:
                    with zf.open(self.annotated_file_path[:-4], mode='r') as f:
                        self.parse_file_object(f)
            else:
                with open(self.annotated_file_path, 'r') as f:
                    self.parse_file_object(f)
        for coding_gene, ch_count in self.ch_counts.items():
            for patient in ch_count.header.patient_columns:
                ch_count["CHCOUNT"] += int(ch_count[patient].value)
        self.save_count_file(count_file_path)

    def parse_file_object(self, file_object):
        encoding = "utf-8"
        first_line = file_object.readline()
        if isinstance(first_line, bytes):
            first_line = first_line.decode(encoding)
        self.header = Header(first_line.strip())
        self.ch_count_header = CodingGeneMutationHeader(self.header, Patients(self.subject_info_file_path))
        self.output_file.write(repr(self.header))
        for line in file_object:
            if isinstance(line, bytes):
                line = line.decode(encoding)
            self.mutation = Mutation(self.header, line.strip(), self.abbreviated_titles)
            if self.mutation["ID"] in self.seen_ids:
                continue
            else:
                self.seen_ids.add(self.mutation["ID"])
            self.mutation.add_frequencies()
            self.mutation.split_info()
            if not self.mutation.rare_variant():
                continue
            self.get_coding_genes()
            for mutation in self.mutations:
                self.parse_new_mutation(mutation)
                self.remove_var_patients(mutation)
                if len(mutation["PATIENTS"]) > 0:
                    self.output_file.write(repr(mutation))

    def get_coding_genes(self):
        """
        Adds the "CODING_GENE" column.
        Each mutation can have multiple coding genes
        This function takes self.mutation, and finds all the coding genes for it
        It the creates copies of self.mutation, each one with one of the coding genes
        """
        self.mutations = []
        coding_genes = self.mutation.get_coding_genes()  # Sets remove duplicates from list
        for coding_gene in coding_genes:
            new_mutation = Mutation(self.mutation)
            new_mutation["CODING_GENE"] = coding_gene
            self.mutations.append(new_mutation)
        self.mutation = None

    def parse_new_mutation(self, mutation):
        """
        The goal of this function is to
        concatenates mutations that affect the same gene
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

        This function modifies the self.ch_counts dictionary.
        If a mutation that affects this coding_gene has already been seen, it modifies that entry with the new
        information.  If not, it creates an entry.
        """
        coding_gene = mutation["CODING_GENE"]
        if coding_gene in self.ch_counts:
            self.ch_counts[coding_gene].parse_new_mutation(mutation)
        else:
            self.ch_counts[coding_gene] = CodingGeneMutation(self.ch_count_header, coding_gene, self.valid_genotypes)
            self.ch_counts[coding_gene].parse_new_mutation(mutation)

    def remove_var_patients(self, mutation):
        """
        This concatenates the patients affected by a mutation
        Beforehand, they were stored as
                patient1    patient2    patient3
        mut1    0|0         1|0         1|1
        Now, they are stored as
                PATIENTS
        mut1    patient2,patient3
        """
        patients_affected_by_mutation = []
        for patient in self.header.patient_columns:
            if mutation.data[patient] != "0|0":
                patients_affected_by_mutation.append(patient)
            del mutation.data[patient]
        mutation["PATIENTS"] = ",".join(patients_affected_by_mutation)

    def save_count_file(self, count_file_path):
        """
        Writes all count files (CodingGeneMutation, CodingGeneMutationHeader) to count_file_path
        """
        with open(count_file_path, 'w') as f:
            f.write(repr(self.ch_count_header))
            for coding_gene, ch_count in self.ch_counts.items():
                f.write(repr(ch_count))

if __name__ == '__main__':
    if len(sys.argv) != 5:
        raise IOError("Command line arguments need to be input_file_path, output_file_path, ch_file_path,"
                      "and subject_info_file_path")
    d = Data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
