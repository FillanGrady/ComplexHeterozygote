#this script makes a new column 'VarID' with an ID that matches the ID created in the corresponding vcf file.

def makeID(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_wID.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    outfile.write('VarID')
    outfile.write('\t')
    outfile.write(line)
    lines = infile.readlines()
    infile.close()
    for line in lines:
        words = line.split('\t')
        varID = words[1] + ':' + words[2] + ',' + words[9] + ',' + words[10]
        outfile.write(varID)
        outfile.write('\t')
        outfile.write(line)
    outfile.close()
            

makeID('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants_clean')





