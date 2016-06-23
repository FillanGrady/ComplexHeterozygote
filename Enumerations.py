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


class GeneEffects(Enum):
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
    START_GAINED = 15
    START_LOST = 16
    STOP_GAINED = 17
    STOP_LOST = 18
    CODON_DELETION = 19
    CODON_INSERTION = 20
    CODON_CHANGE_PLUS_CODON_DELETION = 21
    CODON_CHANGE_PLUS_CODON_INSERTION = 22
    FRAME_SHIFT = 23
    TRANSCRIPT = 24
    SPLICE_SITE_DONOR = 25
    SPLICE_SITE_ACCEPTOR = 26
    INTERGENIC = 27
    INTRAGENIC = 28

    @staticmethod
    def all():
        return [member for name, member in GeneEffects.__members__.items()]

    @staticmethod
    def default():
        return [GeneEffects.CODON_CHANGE_PLUS_CODON_DELETION, GeneEffects.CODON_CHANGE_PLUS_CODON_INSERTION,
                GeneEffects.CODON_DELETION, GeneEffects.CODON_INSERTION, GeneEffects.EXON_DELETED,
                GeneEffects.FRAME_SHIFT, GeneEffects.NON_SYNONYMOUS_CODING, GeneEffects.SPLICE_SITE_ACCEPTOR,
                GeneEffects.SPLICE_SITE_DONOR, GeneEffects.START_LOST, GeneEffects.STOP_GAINED, GeneEffects.START_LOST]

    def __str__(self):
        return self.name
