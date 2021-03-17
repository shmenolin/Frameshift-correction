#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import re
import pandas as pd
from blastx_to_df import *


# In[2]:


#take a sequence and return gaps patterns along with the index of the '-'



class find_gap:
    def __init__(self):
        pass
    def FindGaps(self,string):
        df = pd.DataFrame()
        gap = re.search(r'-', string)
        if gap:
            gaps = re.finditer(r'[ACGTN]{3}-+[ACGTN]{3}', string)
            patterns = []
            positions = []
            for i in gaps:
                pat = i.group()
                ind = i.start()+3
                patterns.append(pat)
                positions.append(ind)
            df['patterns'] = patterns
            df['gap_index'] = positions
            return df

        return df


# In[3]:


def blastn(CDS, db, prot_db, nt_shift, CDS_df):
    old_CDS = CDS
    print ("blastn db is ", db)
    old_start = CDS_df.loc[0, 'start']
    new_start = old_start + nt_shift
    os.system('rm CDS')
    CDS_FILE = open('CDS', 'w') #beacuse -query in blastn requires a file
    CDS_FILE.write(CDS)
    CDS_FILE.close()
    cmd = "blastn -db %s -query CDS -max_target_seqs 1 -perc_identity 80 -outfmt '6 qstart qend qseqid sseqid qseq sseq' > results" %(db)
    os.system(cmd)
    filesize = os.path.getsize("results")
    if filesize > 0:
        df = pd.read_csv("results", sep ='\s+', names = ['qstart', 'qend', 'qseqid', 'sseqid', 'qseq', 'sseq'])
        print("the df is", df, "\n")
    else:
        list1 = [CDS,0]
        return list1

    q_seq = df.loc[0,"qseq"].upper()
    new_q_seq = q_seq
    s_seq = df.loc[0,'sseq'].upper()
    q_start = df.loc[0,"qstart"]-1
    q_end = df.loc[0, "qend"]-1 #make up for py index

    seq = find_gap()
    q_df = seq.FindGaps(q_seq)
    print (q_df)
    list1 = []
    if not q_df.empty:
        for i in range(len(q_df)):
            index = q_df.loc[i, 'gap_index']
            pattern = q_df.loc[i, 'patterns']
            nt_before_gap = re.search(r'[ACGT]-', pattern).group()
            nt_before_gap = nt_before_gap[0]
            nt_after_gap = re.search(r'-[ACGT]', pattern).group()
            nt_after_gap = nt_after_gap[1]
            fix_nt = s_seq[index]
            print ('gap in q_seq:', q_seq[index-3:index+4])
            print ('equivalent s_seq: ', s_seq[index-3:index+4])
            if fix_nt == nt_before_gap or fix_nt == nt_after_gap: #to ensure we have homopolymeric region here, could have done in the above part but doing above does not grab when homopolymer is of length 2nt
                length = pattern.count('-')
                print ("index", index)
                print ("nt_shift", nt_shift)
                print ("old_start", old_start)
                print ("qstart", q_start)
                print ("length", length)
                #print ('equivalent s_seq', s_seq[index-3:index+4])
                pos_new_genome = new_start + q_start+ index
                new_q_seq = new_q_seq[:index] + s_seq[index:index+length] + new_q_seq[index+length:] #replacing the gap with equivalent nt in subject sequence
                change_dict = {'change type': 'insertion', 'position in new genome': pos_new_genome, 'before': '-', 'after': s_seq[index:index+length], 'seq before': q_seq[index-5:index+6], 'seq after': new_q_seq[index-5:index+6]}
                print ('change_dict', change_dict)
                list1.append(change_dict)
    s_df = seq.FindGaps(s_seq)
    print ('gaps in sub seq: ', s_df)


    if not s_df.empty:
        total_change = 0
        for i in range(len(s_df)):
            #if there's an insertion in the query sequence, we expect the qseq to have two similar nts in and before the position
            index = s_df.loc[i, 'gap_index']
            print ('gap in  s_seq: ', s_seq[index-3:index+4])
            print ('equivalent q_seq', q_seq[index-3:index+4])
            if q_seq[index] == q_seq[index-1] or q_seq[index] == q_seq[index+1]:
                #print (i)
                #since we are deleting nts off of the new_q_seq we need to note the shift in nts that occur with deletion
                print ("index", index)
                print ("nt_shift", nt_shift)
                print ("old_start", old_start)
                print ("qstart", q_start)
                #print ("length", length)
                index_new_q_seq = index - total_change
                pattern = s_df.loc[i, 'patterns']
                length = pattern.count('-')
                print ("length", length)
                pos_new_genome = new_start + q_start+ index_new_q_seq
                total_change = total_change + length
                new_q_seq = new_q_seq[:index_new_q_seq] + new_q_seq[index_new_q_seq+length:] #deleting a nt at equivalent position
                change_dict = {'change type': 'deletion', 'position in new genome': pos_new_genome, 'before': q_seq[index:index+length],'after': '', 'seq before': q_seq[index-5:index+6], 'seq after': new_q_seq[index_new_q_seq-5:index_new_q_seq+6]}
                print ('change_dict', change_dict)
                list1.append(change_dict)
    new_q_seq = new_q_seq.replace('-', '') #replacing any remaining '-' with ''
    CDS = CDS[:q_start] + new_q_seq + CDS[q_end+1:] #coz q_end nt is already included in q_seq sequence
    
    if old_CDS == CDS:
        print ("could not be fixed by blastn, cause no gaps in homopolymeric region")
        list = [CDS, 0]
        return list

    print ('CDS from blastn fixing\n', CDS)
    list = [CDS, list1]
    #changed dictionaries are returned as a list of all dictionaries


    #lets return only when it produces single hit against protein db after fixing
    blastx_df = blastx_to_df(CDS, prot_db)
    
    if isinstance(blastx_df, int):
        list = [CDS, 0]
        print ("no hits against protein db after change with blastn")
        return list
    blastx_df = remove_duplicates(blastx_df)
    single_hit = check_single_hit(blastx_df)
    print ('single hit ', single_hit)
    if single_hit == 'yes':
        return list
    else:
        list = [CDS, 0]
        print ("could not be fixed by blastn")
        return list
