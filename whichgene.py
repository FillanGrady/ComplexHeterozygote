#this scripts takes the value for the "effects" column (annotated by SNPeff/SNPsift software) and pulls out the gene name in which the coding mutation is induced.  This inserted into a new column [4], so all columns are shifted over 1.  This may not be necessary - it could be required as part of the user-created input file.  One thing to fix: it pulls out the first instance, so if two genes are affected, only the first one will be pulled out. 

def WHICHGENE(filename):
    infilename = filename +'.txt'
    infile = open(infilename, 'r')
    outfilename = (filename + '_gened.txt')
    outfile = open(outfilename, 'w')
    line = infile.readline()
    words = line.split('\t')
    numwords = len(words)
    for i in range(0, 4):
        outfile.write(words[i])
        outfile.write('\t')
    outfile.write('Coding_gene')
    outfile.write('\t')
    for i in range(4, numwords-1):
        outfile.write(words[i])
        outfile.write('\t')
    outfile.write(words[numwords-1])
    lines = infile.readlines()
    infile.close()
    ibbles = []
    for line in lines:
        words = line.split('\t')
        effects = words[62]
        for i in range(0, 4):
            outfile.write(words[i])
            outfile.write('\t')
        affected = 'unknown'
        effectswords = effects.split(',')
        for effect in effectswords:
            keyword = effect.split('(') 
            if keyword[0] == 'CODON_CHANGE_PLUS_CODON_DELETION':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0] == 'CODON_CHANGE_PLUS_CODON_INSERTION':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0] == 'CODON_INSERTION':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0] == 'CODON_DELETION':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0] == 'NON_SYNONYMOUS_CODING':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'STOP_GAINED':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'STOP_LOST':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'START_GAINED':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'START_LOST':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'FRAME_SHIFT':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'SPLICE_SITE_ACCEPTOR':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'SPLICE_SITE_DONOR':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0]== 'NON_SYNONYMOUS_START':
                genes = effect.split('|')
                affected = genes[5]
                break
            elif keyword[0] == 'EXON_DELETION':
                genes = effect.split('|')
                affected = genes[5]               
        outfile.write(affected)
        outfile.write('\t')
        for i in range(4, numwords-1):
           outfile.write(words[i])
           outfile.write('\t')
        outfile.write(words[numwords-1])
    outfile.close()


WHICHGENE('1kg_chr22_coding_clean_recode1_ACAF_recode2_2pct_noreps')

