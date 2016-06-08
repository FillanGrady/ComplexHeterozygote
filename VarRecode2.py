#in this script, any variant with a MAF > 0.5 is recoded, where 0, 1 and 2 are now 2, 0 and 1, so that the homozygous ref genotype is always 0 and the homozygous var genotype is always 2.
#Col6 [5] is SNP/INDEL
#Col8 [7] is reference allele
#Col9 [8] is var allele
#col [63] is where genotypes begin
#people = 2504

    
def VARrecode2(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_recode2.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    words = line.split('\t')
    for i in range(0, 63):
        outfile.write(words[i])
        outfile.write('\t')
    # AC (allele count), AMax (maximum allele count), and AF (allele frequency) are moved to columns before the genotypes.
    outfile.write('AC')
    outfile.write('\t')
    outfile.write('AMax')
    outfile.write('\t')
    outfile.write('AF')
    outfile.write('\t')
    for i in range(63, 2566):
        outfile.write(words[i])
        outfile.write('\t')
    for i in range(2566, 2567):
        outfile.write(words[i])
        outfile.write('\n')
    lines = infile.readlines()
    infile.close()
    for line in lines:
        words = line.split('\t')
        total = len(words)
        AC = words[2567]
        AMax = words[2568]
        AFnum = words[2569].split('\n')
        AF = AFnum[0]
        for i in range(0,63):
            outfile.write(words[i])
            outfile.write('\t')
        outfile.write(AC)
        outfile.write('\t')
        outfile.write(AMax)
        outfile.write('\t')
        outfile.write(AF)
        outfile.write('\t')
        if float(AF) <= 0.5:
            for i in range(63, 2566):
                if words[i] == 'NA':
                    outfile.write('NA')
                    outfile.write('\t')
                else:
                    outfile.write(words[i])
                    outfile.write('\t')
            for i in range(2566, 2567):
                if words[i] == 'NA':
                    outfile.write('NA')
                    outfile.write('\n')
                else:
                    outfile.write(words[i])
                    outfile.write('\n')           
        else:
            #if AF > 0.5, then alleles should be recoded in the opposite manner
            for i in range(63, 2566):
                if words[i] == 'NA':
                    outfile.write('NA')
                    outfile.write('\t')
                elif words[i] == '0':
                    outfile.write('2')
                    outfile.write('\t')
                elif words[i] == '2':
                    outfile.write('0')
                    outfile.write('\t')
                else:
                    outfile.write(words[i])
                    outfile.write('\t')
            for i in range(2566, 2567):
                if words[i] == 'NA':
                    outfile.write(words[i])
                    outfile.write('\n')
                elif words[i] == '0':
                    outfile.write('2')
                    outfile.write('\n')
                elif words[i] == '2':
                    outfile.write('0')
                    outfile.write('\n')
                else:
                    outfile.write(words[i])
                    outfile.write('\n')
    outfile.close()






VARrecode2('1kg_chr22_coding_clean_recode1_ACAF')









