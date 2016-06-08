#before running this, input file is opened in excel and sorted by genename and then by Subject
#this scripts removes the variants that don't compose possible compound heterozygotes.
#needs to be changed to include homozygous variants too

def cleanVariants(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_clean.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    outfile.write(line)
    lines = infile.readlines()
    infile.close()
    numlines = len(lines)
    written = []
    for i in range(0, numlines):
        written.append('no')
    for i in range(0, numlines):
        words2 = lines[i].split('\t')
        words1 = lines[i-1].split('\t')
        if (words2[0] == words1[0]) and (words2[5] == words1[5]):
            if written[i-1] == 'no':
                outfile.write(lines[i-1])
                written[i-1] = 'yes'
            outfile.write(lines[i])
            written[i] = 'yes'
    outfile.close()
            
        


cleanVariants('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants')
