#Here there are 2 scripts, SepFile and CHCOUNT3v
#File needs to be sorted by col[4] prior to input
#SepFile - this separates the file into information about the variants, and the other file contains the genotype information
#CHCount3V - this takes the output files from SepFile.  Using the "_people" file, for each gene, determines if a person has 2 or more mutations in a gene.
#Then, using info from the "_varinfo" file, variants that compose the compound hets (and all other variants) are output to the "_CHvariants" file and the counts per gene, along with which individual has 2 or more mutations is output in the "_CHcount" file.

def SepFile(filestart):
    filename = (filestart + '.txt')
    outfilename = filestart + '_people.txt'
    outfilename2 = filestart + '_varinfo.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    outfile2 = open(outfilename2, 'w')
    lines= infile.readlines()
    infile.close()
    for line in lines:
        words = line.split('\t')
        numcols = len(words)
        for i in range(0, 66):
            outfile2.write(words[i])
            outfile2.write('\t')
        outfile2.write(words[66])
        outfile2.write('\n')
        outfile.write(words[4])
        outfile.write('\t')
        for i in range(67, numcols-1):
            outfile.write(words[i])
            outfile.write('\t')
        outfile.write(words[numcols-1])


SepFile('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened')


def CHCOUNT3v(filestart, samples):
    filename = filestart + '_people.txt'
    filename2 = filestart + '_varinfo.txt'
    outfilename = filestart + '_CHcount.txt'
    outfilename2 = filestart + '_CHvariants.txt'
    infile = open(filename, 'r')
    infile2 = open(filename2, 'r')
    outfile = open(outfilename, 'w')
    outfile2 = open(outfilename2, 'w')
    line1 = infile.readline()
    line1words = line1.split('\t')
    numwords = len(line1words)
    outfile.write('GENE')
    outfile.write('\t')
    outfile.write('CHCOUNT')
    outfile.write('\t')
    for i in range(1, numwords-1):
        outfile.write(line1words[i])
        outfile.write('\t')
    outfile.write(line1words[numwords-1])
    lines = infile.readlines()
    numlines = len(lines)
    line2_1 = infile2.readline()
    outfile2.write('SUBJECT')
    outfile2.write('\t')
    outfile2.write(line2_1)
    lines2 = infile2.readlines()
    infile.close()
    infile2.close()
    samplecount = []
    samplecountnew = []
    total = 0
    linenum = 0
    for i in range(0, samples):
        samplecount.append(0)
        samplecountnew.append(0)
    genename = 'blank'
    for j in range(0, numlines):
        words = lines[j].split('\t')
        genenamenew = words[0]
        if genenamenew == genename:
            for i in range(1, samples):
                if (words[i] == 'NA'):
                    dummy = 0
                elif ((words[i] == '1') or (words[i] == '2')):
                    samplecount[i-1] = samplecount[i-1] + int(words[i])
                    outfile2.write(line1words[i])
                    outfile2.write('\t')
                    outfile2.write(lines2[j])
                else:
                    samplecountnew[i-1] = samplecountnew[i-1]
        else:
            for i in range(0, samples):
                if samplecount[i] > 1:
                    total = total +1
            outfile.write(genename)
            outfile.write('\t')
            outfile.write(str(total))
            outfile.write('\t')
            for i in range(0, samples):
                outfile.write(str(samplecount[i]))
                outfile.write('\t')
            outfile.write('\n')
            genename = genenamenew
            total = 0
            samplecount = []
            samplecountnew = []
            for i in range(0, samples):
                samplecount.append(0)
                samplecountnew.append(0)
            for i in range(1, samples):
                if (words[i] == 'NA'):
                    dummy = 0
                elif((words[i] == '1') or (words[i] == '2')):
                    samplecountnew[i-1] = samplecountnew[i-1] + int(words[i])
                    outfile2.write(line1words[i])
                    outfile2.write('\t')
                    outfile2.write(lines2[j])
                else:
                    samplecountnew[i-1] = samplecountnew[i-1]
            samplecount = samplecountnew
    for i in range(0, samples):
        if samplecount[i] > 1:
            total = total +1
    outfile.write(genename)
    outfile.write('\t')
    outfile.write(str(total))
    outfile.write('\t')
    for i in range(0, samples):
        outfile.write(str(samplecount[i]))
        outfile.write('\t')
    outfile.close()
    outfile2.close()


CHCOUNT3v('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened', 2504)


