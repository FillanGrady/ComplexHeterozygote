#this scripts removes repeated variants based on dbID (rs#).  I think the "makeID" scripts should be done first and then repeats can be removed based on both ID OR dbID
    
def RemoveRepeats(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_noreps.txt'
    repfilename = filestart + '_reps.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    repfile = open(repfilename, 'w')
    line = infile.readline()
    outfile.write(line)
    lines = infile.readlines()
    numlines = len(lines)
    infile.close()
    for i in range(0, numlines):
        words1 = lines[i-1].split('\t')
        words2 = lines[i].split('\t')
        if (words2[6] == words1[6]):
            repfile.write(lines[i])
        else:
            outfile.write(lines[i])
    outfile.close()



RemoveRepeats('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct')








