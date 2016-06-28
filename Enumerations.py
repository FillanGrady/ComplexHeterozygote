from enum import Enum


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

    @staticmethod
    def from_string(string):
        if "AFR" in string:
            return Ancestry.African
        elif "AMR" in string:
            return Ancestry.American
        elif "EAS" in string:
            return Ancestry.EastAsian
        elif "EUR" in string:
            return Ancestry.European
        elif "SAS" in string:
            return Ancestry.SouthAsian
        else:
            raise KeyError


class Genotypes(Enum):
    """
    Inclusive comparisons;
    OneMutation is for patients with one or more mutations
    """
    CompoundHeterozygotes = 1
    Homozygous = 2
    TwoMutations = 3
    OneMutation = 4


class MutationEffects(Enum):
    """
    What genes are affected by a mutation
    """
    INTRON = 1
    EXON = 2
    EXON_DELETED = 3
    DOWNSTREAM = 4
    UPSTREAM = 5
    SYNONYMOUS_CODING = 6
    SYNONYMOUS_START = 7
    SYNONYMOUS_STOP = 8
    NON_SYNONYMOUS_CODING = 9
    NON_SYNONYMOUS_START = 10
    UTR_3_PRIME = 11
    UTR_5_PRIME = 12
    UTR_3_DELETED = 13
    UTR_5_DELETED = 14
    UTR_3_PREMATURE_START = 15
    UTR_5_PREMATURE_START = 16
    START_GAINED = 17
    START_LOST = 18
    STOP_GAINED = 19
    STOP_LOST = 20
    CODON_DELETION = 21
    CODON_INSERTION = 22
    CODON_CHANGE_PLUS_CODON_DELETION = 23
    CODON_CHANGE_PLUS_CODON_INSERTION = 24
    FRAME_SHIFT = 25
    TRANSCRIPT = 26
    SPLICE_SITE_DONOR = 27
    SPLICE_SITE_ACCEPTOR = 28
    SPLICE_SITE_VARIANT = 29
    INTERGENIC = 30
    INTRAGENIC = 31
    NON_CODING = 32
    MISSENSE = 33
    DISRUPTIVE_INFRAME_DELETION = 34
    DISRUPTIVE_INFRAME_INSERTION = 35

    @staticmethod
    def str_lookup():
        """
        Each annotation software labels the effects of the mutations slightly differently
        This function provides a user-accessible class to allow custom lookups for each mutation_effect
        This returns a dictionary, where the keys are the string representations of a mutation_effect,
        and the values are the mutation_effect

        Ex: d["INTERGENIC_REGION"] = MutationEffects.INTERGENIC
        """
        d = {}
        for name, member in MutationEffects.__members__.items():
            d[str(member)] = member
            d[str(member) + "_VARIANT"] = member  # Frequently, labeled as "SYNONYMOUS_VARIANT" instead of "SYNONYMOUS"
        d["INTERGENIC_REGION"] = MutationEffects.INTERGENIC
        d["DOWNSTREAM_GENE_VARIANT"] = MutationEffects.DOWNSTREAM
        d["UPSTREAM_GENE_VARIANT"] = MutationEffects.UPSTREAM
        d["NON_CODING_EXON_VARIANT"] = MutationEffects.NON_CODING
        d["SPLICE_DONOR_VARIANT"] = MutationEffects.SPLICE_SITE_DONOR
        d["SPLICE_ACCEPTOR_VARIANT"] = MutationEffects.SPLICE_SITE_ACCEPTOR
        d["SPLICE_REGION_VARIANT"] = MutationEffects.SPLICE_SITE_VARIANT
        d["3_PRIME_UTR_VARIANT"] = MutationEffects.UTR_3_PRIME
        d["5_PRIME_UTR_VARIANT"] = MutationEffects.UTR_5_PRIME
        d["SYNONYMOUS_VARIANT"] = MutationEffects.SYNONYMOUS_CODING
        d["FRAME_SHIFT_VARIANT"] = MutationEffects.FRAME_SHIFT
        d["FRAMESHIFT_VARIANT"] = MutationEffects.FRAME_SHIFT
        d["5_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT"] = MutationEffects.UTR_5_PREMATURE_START
        d["3_PRIME_UTR_PREMATURE_START_CODON_GAIN_VARIANT"] = MutationEffects.UTR_3_PREMATURE_START
        return d

    @staticmethod
    def all():
        return [member for name, member in MutationEffects.__members__.items()]

    @staticmethod
    def default():
        return [MutationEffects.CODON_CHANGE_PLUS_CODON_DELETION, MutationEffects.CODON_CHANGE_PLUS_CODON_INSERTION,
                MutationEffects.CODON_DELETION, MutationEffects.CODON_INSERTION, MutationEffects.EXON_DELETED,
                MutationEffects.FRAME_SHIFT, MutationEffects.NON_SYNONYMOUS_CODING, MutationEffects.SPLICE_SITE_ACCEPTOR,
                MutationEffects.SPLICE_SITE_DONOR, MutationEffects.START_LOST, MutationEffects.STOP_GAINED, MutationEffects.START_LOST]

    def __str__(self):
        return self.name
