#this scripts takes the output from VarRecode1 and counts the alleles, so that VarRecode2 can recode genotypes when allele frequency is >0.5
    
def CountAlleles(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_ACAF.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    lineparts = line.split('\n')
    line = lineparts[0]
    #AC (Allele count), AMax (Max allele count) and AF (Allele frequency) are written as new columns at the end.
    outfile.write(line)
    outfile.write('AC')
    outfile.write('\t')
    outfile.write('AMax')
    outfile.write('\t')
    outfile.write('AF')
    outfile.write('\n')
    lines = infile.readlines()
    infile.close()
    for line in lines:
        if line[0] == '#':
            outfile.write(line)
        else:
            words = line.split('\t')
            total = len(words)
            AC = 0
            AF = 0
            AMax = 0
            for i in range(0,total-1):
                outfile.write(words[i])
                outfile.write('\t')
            for i in range(63, total):
                if words[i] == '0':
                    AMax = AMax+2
                elif words[i] == '1':
                    AMax = AMax+2
                    AC = AC + 1
                elif words[i] == '2':
                    AC = AC + 2
                    AMax = AMax+2
                else:
                    AMax = AMax+2
            last = words[total-1].split('\n')
            #AC (Allele count), AMax (Max allele count) and AF (Allele frequency) values are written out.
            outfile.write(last[0])
            outfile.write(str(AC))
            outfile.write('\t')
            AF = AC/AMax
            outfile.write(str(AMax))
            outfile.write('\t')
            outfile.write(str(AF))
            outfile.write('\n')
    outfile.close()





CountAlleles('1kg_chr22_coding_clean_recode1')








