#This script uses a list of variants with genotypes as input
#Col [5] is SNP/INDEL - if a variant is a SNP or an INDEL - this matters because indels are shown as var1:var2 and SNPs are var1var2
#Col [7] is reference allele
#Col [8] is var allele
#col [12] is AF (var)
#need to take reverse alleles into account - because some reference and variant alleles are reversed in the public databases
#col [63] is where genotypes begin
#Function - recodes genotypes as 0,1,and 2 for homozygous reference, heterozygous and homozygous variant, respectively.


    
def VARrecode(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_recode1.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    lines = infile.readlines()
    infile.close()
    #possible genotypes, excluding 'NA'
    hetgenos = ['AC', 'AG', 'AT', 'CG', 'CT', 'CA', 'GC', 'GA', 'GT', 'TA', 'TC', 'TG']
    homogenos = ['AA', 'GG', 'CC', 'TT']
    for line in lines:
        if line[0] == '#':
            outfile.write(line)
        else:
            words = line.split('\t')
            total = len(words)
            type = words[5]
            ref = words[7]
            var = words[8]
            for i in range(0,63):
                outfile.write(words[i])
                outfile.write('\t')
            for i in range(63, total-1):
                if type == 'SNP':
                    if words[i] == ref + ref:
                        outfile.write('0')
                        outfile.write('\t')
                    elif words[i] == ref+var or words[i] == var+ref:
                        outfile.write('1')
                        outfile.write('\t')
                    elif words[i] == var+var:
                        outfile.write('2')
                        outfile.write('\t')
                            #NAs are changed in Recode2
                    elif words[i] == 'NA':
                        outfile.write('NA')
                        outfile.write('\t')
                    elif words[i] in hetgenos:
                        outfile.write('1')
                        outfile.write('\t')
                            #I'm not sure why I did this - although there should be very few cases where this happens, if at all.  Something to figure out and change.
                    elif words[i] in homogenos:
                        outfile.write('2')
                        outfile.write('\t')
                    else:
                        outfile.write(words[i])
                        outfile.write('\t')
                elif type == 'INDEL':
                    if words[i] == ref + ':' + ref:
                        outfile.write('0')
                        outfile.write('\t')
                    elif words[i] == ref + ':' +var or words[i] == var + ':' + ref:
                        outfile.write('1')
                        outfile.write('\t')
                    elif words[i] == var + ':' + var:
                        outfile.write('2')
                        outfile.write('\t')
                    elif words[i] == 'NA':
                        outfile.write('NA')
                        outfile.write('\t')
                    elif ':' in words[i]:
                        half = words[i].split(':')
                        #again, not sure why I assume anything homozygous should be coded as 2.
                        if half[0] == half[1]:
                            outfile.write('2')
                        elif half[0] != half[1]:
                            outfile.write('1')
                        outfile.write('\t')
                    else:
                        outfile.write(words[i])
                        outfile.write('\t')
                else:
                    outfile.write(words[i])
                    outfile.write('\t')
            for i in range(total-1, total):
                wordstart = words[i].split('\n')
                words[i] = wordstart[0]
                if type == 'SNP':
                    if words[i] == ref + ref:
                        outfile.write('0')
                    elif words[i] == ref+var or words[i] == var+ref:
                        outfile.write('1')
                    elif words[i] == var+var:
                        outfile.write('2')
                    elif words[i] == 'NA':
                        outfile.write('NA')
                    elif words[i] in hetgenos:
                        outfile.write('1')
                    elif words[i] in homogenos:
                        outfile.write('2')
                    else:
                        outfile.write(words[i])
                elif type == 'INDEL':
                    if words[i] == ref + ':' + ref:
                        outfile.write('0')
                    elif words[i] == ref + ':' +var or words[i] == var + ':' + ref:
                        outfile.write('1')
                    elif words[i] == var + ':' + var:
                        outfile.write('2')
                    elif words[i] == 'NA':
                        outfile.write('NA')
                    elif ':' in words[i]:
                        half = words[i].split(':')
                        if half[0] == half[1]:
                            outfile.write('2')
                        elif half[0] != half[1]:
                            outfile.write('1')
                    else:
                        outfile.write(words[i])
                else:
                    outfile.write(words[i])
                outfile.write('\n')
                                   
    outfile.close()



VARrecode('1kg_chr22_coding_clean')







