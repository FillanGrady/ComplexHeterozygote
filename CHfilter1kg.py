#before running this, file must be sorted by subject ID and then coding gene
#col[1] = subject ID
#col[2] = gene name
#col[69] = genotype
# for genotype, left side is "dad" allele, right side is "mom" allele

def CHfilter1kg(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_checked.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    temp = open('temp.txt', 'w')
    line = infile.readline()
    outfile.write(line)
    lines = infile.readlines()
    numlines = len(lines)
    infile.close()
    for i in range(0, numlines):
        words1 = lines[i-1].split('\t')
        words2 = lines[i].split('\t')
        if (words2[1] == words1[1]) and (words2[6] == words1[6]):
            temp.write(lines[i-1])
        else:
            if i >= 1:
                temp.write(lines[i-1])
            i = i+2
            temp.close()
            temp = open('temp.txt', 'r')
            genelines = temp.readlines()
            dadyes = 'no'
            momyes = 'no'
            for geneline in genelines:
                words = geneline.split('\t')
                genotype = words[69]
                genotype = genotype.replace('\n', '')
                alleles = genotype.split('|')
                dad = alleles[0]
                mom = alleles[1]
                if dad == '1':
                    dadyes = 'yes'
                if mom == '1':
                    momyes = 'yes'
            if (dadyes == 'yes' and momyes == 'yes'):
                for geneline in genelines:
                    outfile.write(geneline)
            temp.close()
            temp = open('temp.txt', 'w')
    temp.close()
    outfile.close()
            
        



CHfilter1kg('1kg_chr19_CHvariants_phased')
CHfilter1kg('1kg_chr20_CHvariants_phased')
CHfilter1kg('1kg_chr21_CHvariants_phased')
CHfilter1kg('1kg_chr22_CHvariants_phased')

