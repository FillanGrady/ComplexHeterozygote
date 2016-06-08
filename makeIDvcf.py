#This scripts takes a vcf file and makes an ID that is can be matched with the IDs in the CH output file. Where VarID is composed of 'Chr:bppos,refallele,varallele'
#Before this, vcfs are opened and all headers are removed.

def makeID(filestart):
    filename = filestart + '.vcf'
    outfilename = filestart + '_wID.vcf'
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
        leftpos = int(words[1]) - 1
        varID = 'chr' + words[0] + ':' + str(leftpos) + ',' + words[3] + ',' + words[4]
        outfile.write(varID)
        outfile.write('\t')
        outfile.write(line)
    outfile.close()
            



makeID('chr22.ANNOTATED.db.eff_coding')







