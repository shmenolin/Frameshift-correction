#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import sys
import re
import pandas as pd
import argparse
#from blastn import *
from blastx_to_df import *
from extract_CDS import *
from iterative_frmeshift_fix2 import blastx_CDS
#import premature_stop
import csv
from blastn import blastn
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('max_colwidth', 50)

# In[2]:

def get_ntShift(CDS_df, full_CDS, nt_shift):
    start = CDS_df.loc[0,"start"] #index in genome is py index + 1 but already accounted for
    stop = CDS_df.loc[0,"stop"]#this is not accounted for py index yet
    old_CDS_length = stop- start #stop is not fixed to py index
    new_CDS_length = len(full_CDS)
    nt_shift = nt_shift + new_CDS_length - old_CDS_length
    print("the fixed CDS increased by: ",new_CDS_length - old_CDS_length, "and total shift is: ",  nt_shift, "\n")
    return nt_shift

def integrate(CDS_df, full_CDS, new_genome,nt_shift, genome):
    start = CDS_df.loc[0,"start"] #index in genome is py index + 1 but already accounted for
    new_start = int(start) +  nt_shift
    stop = CDS_df.loc[0,"stop"]#this is not accounted for py index yet
    new_stop = int(stop) + nt_shift
    integrated = False
    #if new_genome[new_start:new_start+5] == full_CDS[:5]:
    print ("the region from new start to new stop in new genome is: (should be same as new CDS below) as this region is replaced by new CDS ", new_genome[new_start:new_stop])
    print ("the new full CDS is \n", full_CDS)
    print ('integrating the CDS to genome')
    new_genome = new_genome[:new_start] + full_CDS + new_genome[new_stop:]
    integrated = True
    #print ("the region from new start to new stop in new genome is: (should be same as new CDS below) as this region is replaced by new CDS ", new_genome[new_start:new_stop])
    #print ("the new full CDS is \n", full_CDS)
    print ("earlier it was +- 3 nts from CDS : \n", genome[start-3:stop+4], "\n")
    print ("after being fixed +- 3 nts from CDS: \n", new_genome[new_start-3:new_stop+4])
    #print ('not integrating. Looks like a change to the same region!')
    list = [new_genome, integrated]
    return list


# In[ ]:

def split_string(string):
    string2 = ''
    for i in range(int(len(string)/70)+1):
        if i == int(len(string)/70):
            string2 = string2 + string[i*70:]
            return string2

        string2 = string2 + string[i*70:(i+1)*70] + "\n"


def main():
    #first change the truncated protein list to a dataframe
    parser = argparse.ArgumentParser(description= "Correct frameshifting indels in genes in a pacbio sequencing assembled genome")

    parser.add_argument("Truncated_proteins_file", type = str, help = 'csv file containing the "query_num", "qlen", "slen","protein_name", "start", "stop", "slen-qlen", "qlen_covered" of each CDS')

    parser.add_argument("genome_file", type = str, help = 'genome file in fasta format in which indel has to be fixed')

    parser.add_argument("protein_database", type = str, help = 'name of protein database which will be used for blastx of the CDS sequences')
    parser.add_argument("nt_database",type = str, help = 'name of nucleotide database which will be used for blastn of the CDS sequences')
    #parser.add_argument("--premature", help = "pass the option if you want to look for premature stop codons as well")
    args = parser.parse_args()
    # if args.premature:
    #     fix_premature = 'yes'
    # else:
    #     fix_premature = 'No'

    file = args.Truncated_proteins_file
    db = args.protein_database
    nt_db = args.nt_database
    genome_file_name = args.genome_file
    #file = '/Users/menolinsharma/Documents/correct_frameshifts/indel_correction/scripts/prokka_files/reset_startK03D_contig1.fa/sorted_truncated_proteins.csv' #
    #db = '/Users/menolinsharma/Documents/correct_frameshifts/indel_correction/db/AllBradyProteins_ncbi.faa' #args.protein_database
    #nt_db ='/Users/menolinsharma/Documents/correct_frameshifts/indel_correction/db/all_cds_brady.fa' #args.nt_database
    #genome_file_name = '/Users/menolinsharma/Documents/correct_frameshifts/indel_correction/scripts/prokka_files/reset_startK03D_contig1.fa/eset_sta.fna' #args.genome_file
    genome_file1 = open(genome_file_name, 'r')
    genome_file = genome_file1.readlines()
    genome_file1.close()
    new_genome_file = open("new_iterative_trial1", 'w+')
    CDSes_fixed = open('fixed_CDS_files.fa', 'w+')
    header = ["query_num", 'qlen', 'slen', 'protein_name', 'start', 'stop', "slen_Minus_qlen", "qlen_covered"]
    truncated_proteins_df = pd.read_csv(file, names = header) #list of proteins should be sorted from position 1 to last in the genome
    nt_shift = 0
    change_info_list = []
    genome = ''
    for lines in genome_file:
        if not re.match(r'>', lines):
            genome = genome+ lines.strip() #making big genome string

    new_genome = genome #changes will occur in new genome
    not_corrected = []
    single_hit_afterExtension = []
    no_hits_db = []
    fixed_proteins = []
    proteinsWithStopCodons = []

    CDS_df = extract_CDS(truncated_proteins_df, genome)
            #now we have a dataframe containing the names of all CDS, their sequence and their start and stop positions in genome
            #The adjacent CDSes are joined together


    print("\n the CDS df is\n" + str(CDS_df))

    for i in range(len(CDS_df)):
        #print (CDS_df.index, i)
        #print ('\n CDS_df \n',CDS_df.loc[[i], :]) #i is given in list so that it returns a dataframe instead of series
        CDS_old_genome = str(CDS_df.loc[i,'CDS'])
        CDS =new_genome[CDS_df.loc[i,'start']+nt_shift:CDS_df.loc[i,'stop']+nt_shift]

        print(("\n\n-----------------------------------------\nprocessing protein name:  " + str(CDS_df.loc[i,"protein_name"])))

        print("\nCDS before being extended or fixed\n" + CDS + "\nlength of CDS before being fixed:  "+ str(len(CDS))+"\n")
        print ('CDS from old genome', CDS_old_genome)
        print ('CDS from new genome', CDS)

        #print ("CDS", CDS, "\n")
        df_CDS = CDS_df.loc[[i], :]
        df_CDS.index = [0] #set the index to 0
        print("CDS dataframe for the entry\n" , str(df_CDS))
        print(("\nfirst blastx_df:" ))
        blastx_df = blastx_to_df(CDS, db)
        ProteinName = CDS_df.loc[i,"protein_name"]
        if isinstance(blastx_df, int):
            print("Did not produce any hits against protein db, moving to next entry \n")
            no_hits_db.append([ProteinName, CDS])
            continue

        blastx_df = remove_duplicates(blastx_df)
        single_hit = check_single_hit(blastx_df)
        if single_hit == 'yes':
            print("\nsingle hit against protein db received, will try to extend the region if the gene was found to be truncated\n")
            extended_df_CDS = extend_truncated(blastx_df, CDS_df.loc[[i], :], genome, i)
            CDS1 = extended_df_CDS.loc[0,"CDS"]

            if len(CDS1) == len(CDS):
                print("\ncould not be extended as it seems not to be truncated\n")
                continue

            CDS = CDS1
            CDS=new_genome[extended_df_CDS.loc[0,"start"]+nt_shift:extended_df_CDS.loc[0,"stop"]+nt_shift] #new_change_here
            print("blastx result for extended CDS \n")
            blastx_df = blastx_to_df(CDS, db )
            if isinstance(blastx_df, int):
                print("Did not produce any hits against protein db, moving to next entry \n")
                continue

            blastx_df = remove_duplicates(blastx_df)
            single_hit = check_single_hit(blastx_df)

            if single_hit == 'yes': #if it gives a single hit even after extension, move to another CDS. Leave this one
                print("it produced only single hit even after extending will try to find premature stop codons now if --premature option is provided")
                single_hit_afterExtension.append(extended_df_CDS)
                # if fix_premature == 'yes':
                #     new_genome, changes_dict_list = premature_stop.main(extended_df_CDS, blastx_df, nt_shift, CDS, new_genome, genome)
                #     if changes_dict_list != 0:
                #         for i in changes_dict_list:
                #             change_info_list.append(i)
                #         proteinsWithStopCodons.append(df_CDS)



            else:
                print("more than 1 hits found against the protein db. Trying to fix the region using blastn against the database now... \n")
                changed_list = blastn(CDS, nt_db, db, nt_shift, extended_df_CDS)
                change_dicts = changed_list[1]
                if not isinstance(change_dicts, int): #it is 0 if no correction has been made
                    print ('blastn is approved!')
                    fixed_CDS = changed_list[0]
                    fixed_proteins.append(ProteinName)
                    new_genome, integrated = integrate(extended_df_CDS, fixed_CDS, new_genome, nt_shift, genome)
                    if integrated == True:
                        nt_shift = get_ntShift(extended_df_CDS, fixed_CDS, nt_shift)
                        CDSes_fixed.write('>'+ProteinName+'\n'+fixed_CDS+'/n')
                        for changed_dict in change_dicts:
                            change_info_list.append(changed_dict)
                    continue



                print("more than 1 hits found against the protein db. Trying to fix the region iteratively now using blastx... \n")
                #input should be CDS file
                #output should be changed CDS sequence along with the positions of change and change types
                #output should be zero if it still produced multiple hits in the iterative script at the last check point
                header_CDS = '>'+ CDS_df.loc[i,"protein_name"] + '\n'
                CDS = header_CDS + CDS
                CDS_file = open('CDSfile', 'w+')
                CDS_file.write(CDS)
                CDS_file.close()
                CDS_start = extended_df_CDS.loc[0, 'start']
                if len(blastx_df) > 2:
                    print ("more than two hits not fixed at the moment")
                    not_corrected.append(ProteinName)
                    continue
                changed_list = blastx_CDS('CDSfile', CDSes_fixed, blastx_df, CDS_start, nt_shift, db, nt_db)
                if isinstance(changed_list, int):
                    not_corrected.append(ProteinName)
                    continue
                else:
                    fixed_CDS = changed_list[0]
                    change_dict = changed_list[1]
                    fixed_proteins.append(ProteinName)
                    #change_info_list.append(change_dict)
                    new_genome, integrated = integrate(extended_df_CDS, fixed_CDS, new_genome, nt_shift, genome)
                    if integrated == True:
                        nt_shift = get_ntShift(extended_df_CDS, fixed_CDS, nt_shift)
                        change_info_list.append(change_dict)



        else:
            CDS=new_genome[df_CDS.loc[0,"start"]+nt_shift:df_CDS.loc[0,"stop"]+nt_shift]
            print("there were multiple hits to protein db, first we will try blastn and see if it fixes the indel as it will insert/delete in correct position\n")
            #print("more than 1 hits found against the protein db. Trying to fix the region using blastn against the database now... \n")
            changed_list = blastn(CDS, nt_db, db, nt_shift, df_CDS)
            change_dicts = changed_list[1]
            if not isinstance(change_dicts, int):
                print ('blastn is approved!')
                fixed_CDS = changed_list[0]
                fixed_proteins.append(ProteinName)
                new_genome, integrated = integrate(df_CDS, fixed_CDS, new_genome, nt_shift, genome)
                if integrated == True:
                    nt_shift = get_ntShift(df_CDS, fixed_CDS, nt_shift)
                    CDSes_fixed.write('>'+ProteinName+'\n'+fixed_CDS+'/n')
                    for changed_dict in change_dicts:
                        change_info_list.append(changed_dict)
                continue

            header_CDS = '>'+ CDS_df.loc[i,"protein_name"] + '\n'
            CDS = header_CDS + CDS
            CDS_file = open('CDSfile', 'w+')
            CDS_file.write(CDS)
            CDS_file.close()
            CDS_start = df_CDS.loc[0, 'start']
            if len(blastx_df) > 2:
                print ("more than two hits not fixed at the moment")
                not_corrected.append(ProteinName)
                continue
            changed_list = blastx_CDS('CDSfile', CDSes_fixed, blastx_df, CDS_start, nt_shift, db, nt_db)
            if isinstance(changed_list, int):
                not_corrected.append(ProteinName)
                continue
            else:
                fixed_CDS = changed_list[0]
                change_dict = changed_list[1]
                fixed_proteins.append(ProteinName)
                #change_info_list.append(change_dict)
                new_genome, integrated = integrate(df_CDS, fixed_CDS, new_genome, nt_shift, genome)
                if integrated == True:
                    nt_shift = get_ntShift(df_CDS, fixed_CDS, nt_shift)
                    change_info_list.append(change_dict)
    
    
    new_genome_file.write(">frameshift_fixed_genome_"+ genome_file_name[:7]+ "\n")
    new_genome_70char = split_string(new_genome)
    new_genome_file.write(new_genome_70char + "\n")
    new_genome_file.close
    print ("the change infos are: \n", change_info_list)
    if len(change_info_list) >=1:
        for i in change_info_list:
            print (i)
    keys = change_info_list[0].keys()
    with open("changed_positions.csv", "w+") as csv_file:
            dict_writer = csv.DictWriter(csv_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(change_info_list)
    csv_file.close()


    print ("\n---------------------------\nthe CDSes that gave no hits against protein database")
    for i in no_hits_db:
        print (i, "\n")
    print ("\n---------------------------\nproteins with stop codons are")

    for i in proteinsWithStopCodons:
        print (i)


    print ("\n---------------------------\nthe CDSes giving single hit after extension as well")
    for CDSdf in single_hit_afterExtension:
        print (CDSdf, "\n")


    print ("\n---------------------------\nthe blastx results for CDS that gave multiple hits and could not be corrected")
    for blastxdf in not_corrected:
        print (blastxdf, "\n")


    print("\n---------------------------\nthe corrected proteins with name and CDS")

    for i in fixed_proteins:
        print (i, "\n")

    CDSes_fixed.close()
    #new_genome_file.write(">frameshift_fixed_genome_"+ genome_file_name[:7]+ "\n")
    #new_genome_70char = split_string(new_genome)
    #new_genome_file.write(new_genome_70char + "\n")
    #new_genome_file.close



# In[ ]:


if __name__=="__main__":
    main()
