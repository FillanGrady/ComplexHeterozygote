import time
import datetime
import os
from MultiProcessing import MultiProcessing, conventional
"""
This script uses a list of variants with genotypes as input
Col [5] is SNP/INDEL - if a variant is a SNP or an INDEL - this matters because indels are shown as var1:var2 and SNPs are var1var2
Col [7] is reference allele
Col [8] is var allele
col [12] is AF (var)
need to take reverse alleles into account - because some reference and variant alleles are reversed in the public databases
col [63] is where genotypes begin
Function - recodes genotypes as 0,1,and 2 for homozygous reference, heterozygous and homozygous variant, respectively.
"""


def var_recode(line):
    hetgenos = ['AC', 'AG', 'AT', 'CG', 'CT', 'CA', 'GC', 'GA', 'GT', 'TA', 'TC', 'TG']
    homogenos = ['AA', 'GG', 'CC', 'TT']
    output_line = ""
    if line[0] == '#':
        output_line += line
    else:
        words = line.split('\t')
        total = len(words)
        type = words[5]
        ref = words[7]
        var = words[8]
        for i in range(0, 63):
            output_line += words[i] + '\t'
        for i in range(63, total - 1):
            if type == 'SNP':
                if words[i] == ref + ref:
                    output_line += '0\t'
                elif words[i] == ref + var or words[i] == var + ref:
                    output_line += '1\t'
                elif words[i] == var + var:
                    output_line += '2\t'
                    # NAs are changed in Recode2
                elif words[i] == 'NA':
                    output_line += 'NA\t'
                elif words[i] in hetgenos:
                    output_line += '1\t'
                    # I'm not sure why I did this - although there should be very few cases where this happens, if at all.  Something to figure out and change.
                elif words[i] in homogenos:
                    output_line += '2\t'
                else:
                    output_line += words[i] + '\t'
            elif type == 'INDEL':
                if words[i] == ref + ':' + ref:
                    output_line += '0\t'
                elif words[i] == ref + ':' + var or words[i] == var + ':' + ref:
                    output_line += '1\t'
                elif words[i] == var + ':' + var:
                    output_line += '2\t'
                elif words[i] == 'NA':
                    output_line += 'NA\t'
                elif ':' in words[i]:
                    half = words[i].split(':')
                    # again, not sure why I assume anything homozygous should be coded as 2.
                    if half[0] == half[1]:
                        output_line += '2\t'
                    else:
                        output_line += '1\t'
                else:
                    output_line += words[i] + '\t'
            else:
                output_line += words[i] + '\t'
        for i in range(total - 1, total):
            wordstart = words[i].split('\n')
            words[i] = wordstart[0]
            if type == 'SNP':
                if words[i] == ref + ref:
                    output_line += '0'
                elif words[i] == ref + var or words[i] == var + ref:
                    output_line += '1'
                elif words[i] == var + var:
                    output_line += '2'
                elif words[i] == 'NA':
                    output_line += 'NA'
                elif words[i] in hetgenos:
                    output_line += '1'
                elif words[i] in homogenos:
                    output_line += '2'
                else:
                    output_line += words[i]
            elif type == 'INDEL':
                if words[i] == ref + ':' + ref:
                    output_line += '0'
                elif words[i] == ref + ':' + var or words[i] == var + ':' + ref:
                    output_line += '1'
                elif words[i] == var + ':' + var:
                    output_line += '2'
                elif words[i] == 'NA':
                    output_line += 'NA'
                elif ':' in words[i]:
                    half = words[i].split(':')
                    if half[0] == half[1]:
                        output_line += '2'
                    elif half[0] != half[1]:
                        output_line += '1'
                else:
                    output_line += words[i]
            else:
                output_line += words[i]
    return output_line

if __name__ == '__main__':
    start_time = time.time()
    conventional(input_file_path='AnnotatedTest.txt', output_file_path='fillan.txt', f=var_recode)
    print "Time to execute: %s" % datetime.timedelta(seconds=time.time() - start_time)
