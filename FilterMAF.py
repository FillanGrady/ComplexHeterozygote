#this script removes variants with minor allele frequency > 2% in any population - this is something that should be specified by the user - what MAF cutoff is wanted and based on which population MAF
#Col [13] is 1kg GMAF
#Col [65] is Allele frequency in this set of genomes
#col[14] is AFR 1kg MAF
#col[15] is AMR 1kg MAF
#col[21] is EAS 1kg MAF
#col[23] is EUR 1kg MAF
#col[34] is EUR 1kg MAF
#col[39] is EVS MAFs, comma separated - these are in actual percentages, so 2% MAF is 2.0 not 0.02
#mysterious hidden quotes in EVS MAFs - they are removed in this code.  Only in some chromosomes

    
def FilterMAF(filestart):
    filename = filestart + '.txt'
    outfilename = filestart + '_2pct.txt'
    infile = open(filename, 'r')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    outfile.write(line)
    lines = infile.readlines()
    infile.close()
    for line in lines:
        EVSAFs = [3,3,3]
        words = line.split('\t')
        EVSAFs = words[39].split(',')
        #if MAF = 0, it is annotated as '-'
        if EVSAFs[0] == '-':
            EUEVS = 0
            AAEVS = 0
        else:
            EVSAFs[0] = EVSAFs[0].replace('"','')
            EUEVS = float(EVSAFs[0])
            EVSAFs[1] = EVSAFs[1].replace('"','')
            AAEVS = float(EVSAFs[1])
        EA1KG = words[23]
        if EA1KG == '-':
            EA1KGAF = 0
        else:
            EA1KGAF = float(EA1KG)
        AFR1KG = words[14]
        if AFR1KG == '-':
            AFR1KGAF = 0
        else:
            AFR1KGAF = float(AFR1KG)
        AF= float(words[65])
        GMAF = float(words[13])
        #this is something that I would like to depend on specification by user - some may only want to consider AFs from specific populations
        #but in this case, I am removing all variants with MAF > 2% in any population
        if (((AF <= 0.02) or (AF >= 0.98)) and (GMAF <= 0.02) and (EUEVS <= 2.0) and (AAEVS <= 2.0) and (EA1KGAF <= 0.02 or EA1KGAF >= 0.98) and (AFR1KGAF <= 0.02 or AFR1KGAF >= 0.98)):
            outfile.write(line)

    outfile.close()



FilterMAF('1kg_chr22_coding_clean_recode1_ACAF_recode2')





