#in CH variants file:
#variant ID = column [0]
#subject = column[1]
#dbID = col[9]
#in vcf:dbID = col[3]
#infile1 = CHvariants
#infile2 = ColumnIndex (for vcf file)
#infile3 = vcf file
#before running this, variant files need to be sorted by left flank bp pos.
#modified so that with each successive lookup, starts at last row where positive match was found - that way, don't start at the beginning of the file everytime.


def CHlookup_1kg(filename1, filename2, filename3, filename4):
    infile1 = open(filename1, 'r')
    outfile = open(filename4, 'w')
    line1 = infile1.readline()
    words = line1.split('\t')
    numwords = len(words)
    for i in range(0, numwords-1):
        outfile.write(words[i])
        outfile.write('\t')
    lastword = words[numwords-1].split('\n')
    outfile.write(lastword[0])
    outfile.write('\t')
    outfile.write('vcf_genotype')
    outfile.write('\n')
    lines1 = infile1.readlines()
    infile1.close()
    rowstart = 0
    for line in lines1:
        words = line.split('\t')
        for i in range(0, numwords-1):
            outfile.write(words[i])
            outfile.write('\t')
        lastword = words[numwords-1].split('\n')
        outfile.write(lastword[0])
        outfile.write('\t')
        varID = words[0]
        varID = varID.replace('"','')
        dbID = words[9]
        dbID = dbID.replace('"', '')
        subject = words[1].split('.NA')
        subject = subject[0]
        infile2 = open(filename2, 'r')
        lines2 = infile2.readlines()
        infile2.close()
        for line2 in lines2:
            words2 = line2.split('\t')
            if subject == words2[0]:
                colnum = int(words2[1])
                infile3 = open(filename3, 'r')
                lines3 = infile3.readlines()
                numlines3 = len(lines3)
                for i in range (rowstart, numlines3):
                    words3 = lines3[i].split('\t')
                    if dbID == words3[3]:
                        found = 'TRUE'
                        outfile.write(words3[colnum])
                        rowstart = i
                        outfile.write('\n')
                        infile3.close()
                        break
                    elif varID == words3[0]:
                        found = 'TRUE'
                        outfile.write(words3[colnum])
                        rowstart = i
                        outfile.write('\n')
                        infile3.close()
                        break
                    else:
                        found = 'FALSE'
                if found == 'FALSE':
                    outfile.write('\n')
                break
            else:
                dummy = 0
    outfile.close()



CHlookup_1kg('1kg_chr14_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants_clean_wID.txt', 'Column_Index.txt', 'chr14.ANNOTATED.db.eff_coding_wID.vcf', '1kg_chr14_CHvariants_phased.txt')
CHlookup_1kg('1kg_chr15_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants_clean_wID.txt', 'Column_Index.txt', 'chr15.ANNOTATED.db.eff_coding_wID.vcf', '1kg_chr15_CHvariants_phased.txt')
CHlookup_1kg('1kg_chr16_coding_clean_recode1_ACAF_recode2_2pct_noreps_gened_CHvariants_clean_wID.txt', 'Column_Index.txt', 'chr16.ANNOTATED.db.eff_coding_wID.vcf', '1kg_chr16_CHvariants_phased.txt')


                
                
                
                
    
    
    
