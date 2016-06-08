#before running this, file must be sorted by coding gene and then subject ID.

def TrueCHcount(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_counts.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    lines = infile.readlines()
    infile.close()
    count = 0
    genename = 'dummy'
    subject = 'dummy'
    for line in lines:
        words = line.split('\t')
        subjectnew = words[1]
        genenamenew = words[7]
        if genenamenew == genename:
            if subjectnew != subject:
                count = count + 1
                subject = subjectnew
        elif genenamenew != genename:
            outfile.write(genename)
            outfile.write('\t')
            outfile.write(str(count))
            outfile.write('\n')
            genename = genenamenew
            subject = subjectnew
            count = 1

    outfile.close()

TrueCHcount('Epi4k_Chr1to22_CHgenes_April2015')
