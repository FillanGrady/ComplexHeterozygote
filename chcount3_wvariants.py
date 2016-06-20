import re
#Here there are 2 scripts, SepFile and CHCOUNT3v
#File needs to be sorted by col[4] prior to input
#SepFile - this separates the file into information about the variants, and the other file contains the genotype information
#CHCount3V - this takes the output files from SepFile.  Using the "_people" file, for each gene, determines if a person has 2 or more mutations in a gene.
#Then, using info from the "_varinfo" file, variants that compose the compound hets (and all other variants) are output to the "_CHvariants" file and the counts per gene, along with which individual has 2 or more mutations is output in the "_CHcount" file.

def SepFile(filestart):
    """
    This function splits the input file into two files
    Columns 4, 67- (Coding gene, Patient numbers) go into _varinfo
    Columns 1-3, 5-66 go into _people.txt
    """
    filename = (filestart + '.txt')
    outfilename = filestart + '_people.txt'
    outfilename2 = filestart + '_varinfo.txt'
    outfile = open(outfilename, 'w')
    outfile2 = open(outfilename2, 'w')
    with open(filename, 'r') as f:
        for line in re.split('\r|\n', f.read()):
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
            outfile.write('\n')


SepFile('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened')


def CHCOUNT3v(filestart, samples):
    filename = filestart + '_people.txt'
    filename2 = filestart + '_varinfo.txt'
    outfilename = filestart + '_CHcount.txt'
    outfilename2 = filestart + '_CHvariants.txt'

    input_people_file = open(filename, 'r')
    people_lines = re.split('\r|\n', input_people_file.read().strip())
    input_people_file.close()

    input_var_file = open(filename2, 'r')
    var_lines = re.split('\r|\n', input_var_file.read().strip())
    input_var_file.close()

    count_file = open(outfilename, 'w')
    var_file = open(outfilename2, 'w')
    line1words = people_lines[0].split('\t')
    numwords = len(line1words)
    count_file.write('GENE')
    count_file.write('\t')
    count_file.write('CHCOUNT')
    count_file.write('\t')
    for i in range(1, numwords-1):
        count_file.write(line1words[i])
        count_file.write('\t')
    count_file.write(line1words[numwords-1])
    var_file.write('SUBJECT')
    var_file.write('\t')
    var_file.write(var_lines[0])
    samplecount = []
    samplecountnew = []
    total = 0
    for i in range(0, samples):
        samplecount.append(0)
        samplecountnew.append(0)
    genename = 'blank'
    for j in range(1, len(people_lines)):  # Ignore header line
        words = people_lines[j].split('\t')
        genenamenew = words[0]
        if genenamenew == genename:  # If this line's gene name is the same as the last one's
            print len(words)
            for i in range(1, samples):
                if words[i] == '1' or words[i] == '2':
                    samplecount[i-1] += int(words[i])
                    var_file.write(line1words[i])
                    var_file.write('\t')
                    var_file.write(var_lines[j])
        else:
            for i in range(0, samples):
                if samplecount[i] > 1:
                    total += 1
            count_file.write(genename)
            count_file.write('\t')
            count_file.write(str(total))
            count_file.write('\t')
            for i in range(0, samples):
                count_file.write(str(samplecount[i]))
                count_file.write('\t')
            count_file.write('\n')
            genename = genenamenew
            total = 0
            samplecount = []
            samplecountnew = []
            for i in range(0, samples):
                samplecount.append(0)
                samplecountnew.append(0)
            for i in range(1, samples):
                if words[i] == '1' or words[i] == '2':
                    samplecountnew[i-1] += int(words[i])
                    var_file.write(line1words[i])
                    var_file.write('\t')
                    var_file.write(var_lines[j])
            samplecount = samplecountnew
        var_file.write('\n')
    for i in range(0, samples):
        if samplecount[i] > 1:
            total += 1
    count_file.write(genename)
    count_file.write('\t')
    count_file.write(str(total))
    count_file.write('\t')
    for i in range(0, samples):
        count_file.write(str(samplecount[i]))
        count_file.write('\t')
    count_file.close()
    var_file.close()


CHCOUNT3v('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened', 2504)


