#this scripts matches up IDs between the CHvariants list and the original vcf file - the output contains the allele information (which chromosome) for each variant.
#this takes too long - should be changed to possibly where the query starts in the middle and moves up or down..
#also, I think that the overall program (compilation of scripts) could just start with the annotated vcf.  

def MATCHUP(filename1, filename2, filename3):
    infile1 = open(filename1, 'r')
    infile2 = open(filename2, 'r')
    outfile = open(filename3, 'w')
    line1 = infile1.readline()
    lines1 = infile1.readlines()
    infile1.close()
    line2 = infile2.readline()
    outfile.write(line2)
    lines2 = infile2.readlines()
    infile2.close()
    for line in lines1:
        match = 'no'
        words = line.split('\t')
        varID = words[0]
        varID = varID.replace('"','')
        dbID = words[9]
        dbID = dbID.replace('"', '')
        for line in lines2:
            words2 = line.split('\t')
            varID2 = words2[0]
            varID2 = varID2.replace('"','')
            dbID2 = words2[3]
            dbID2 = dbID2.replace('"','')
            if (dbID == dbID2):
                outfile.write(line)
                match = 'yes'
            elif (varID == varID2):
                outfile.write(line)
                match = 'yes'
            if match == 'yes':
                break
    outfile.close()


MATCHUP('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants_clean_wID.txt', 'chr22.ANNOTATED.db.eff_coding_wID.vcf', '1kg_chr22matchup.txt')
