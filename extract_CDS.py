#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import re
import pandas as pd


# In[2]:


#a function to concatenate start and stop positions of a CDS until all the adjacent ones are joined together
#should result in a dataframe containind protein name, CDS, start stop positions in genome

def extract_CDS(truncatedProteins_df, genome):
    protein_names1 = truncatedProteins_df.loc[:,'protein_name'] #at first lets fix for only single protein split into two
    print ("protein names \n", protein_names1)
    protein_names=[]
    for i in range(len(protein_names1)):
        name = re.search(r'.+', protein_names1[i]).group()
        protein_names.append(name)

    CDS_df = pd.DataFrame(columns=['protein_name', 'CDS', 'start', 'stop'])
    for i in range(len(protein_names)):
        if i == 0:
            protein = protein_names[i]
            start =int(truncatedProteins_df.loc[i, 'start']) -1
            stop = int(truncatedProteins_df.loc[i,'stop'])
            if len(protein_names) == 1:
                CDS = genome[start:stop]
                new_data = {'protein_name': protein, 'CDS':CDS, 'start': start, 'stop': stop}
                CDS_df = CDS_df.append(new_data, ignore_index = True)

        elif i==len(protein_names) - 1:
            if protein_names[i] == protein_names[i-1]:
                stop = int(truncatedProteins_df.loc[i,'stop'])
                CDS = genome[start:stop]
                new_data = {'protein_name': protein, 'CDS':CDS, 'start': start, 'stop': stop}
                CDS_df = CDS_df.append(new_data, ignore_index = True)
            else:
                CDS = genome[start:stop]
                new_data = {'protein_name': protein, 'CDS':CDS, 'start': start, 'stop': stop}
                CDS_df = CDS_df.append(new_data, ignore_index = True)
                start = int(truncatedProteins_df.loc[i, 'start']) -1
                stop = int(truncatedProteins_df.loc[i,'stop'])
                CDS = genome[start:stop]
                protein = protein_names[i]
                new_data = {'protein_name': protein, 'CDS':CDS, 'start': start, 'stop': stop}
                CDS_df = CDS_df.append(new_data, ignore_index = True)



        elif protein_names[i] != protein_names[i-1]:
            CDS = genome[start:stop]
            new_data = {'protein_name': protein, 'CDS':CDS, 'start': start, 'stop': stop}
            CDS_df = CDS_df.append(new_data, ignore_index = True)
            protein = protein_names[i]
            start = int(truncatedProteins_df.loc[i, 'start']) -1
            stop = int(truncatedProteins_df.loc[i,'stop'])

        elif protein_names[i] == protein_names[i-1]:
            stop = int(truncatedProteins_df.loc[i,'stop'])

    print ('\nAll truncated proteins dataframe\n', CDS_df)
    return CDS_df


# In[ ]:

#to see if truncated protein could use some help with extension
#done only when single hit function produces yes

def extend_truncated(blastx_df, CDS_df, genome, i):
    #we suppose this is giving us single hit
    start = CDS_df.loc[i, 'start'] #is already -1 and made up for py index
    stop = CDS_df.loc[i,'stop'] #no -1 coz python excludes the last character anyway
    #CDS = CDS_df.loc[i,'CDS']
    #print ("length of CDS", len(CDS),"\n")
    name = CDS_df.loc[i,'protein_name']
    slen = blastx_df.loc[0,'slen']
    sstart = blastx_df.loc[0,'sstart']
    send = blastx_df.loc[0,'send']

    if sstart < send and ((sstart - 0) > 10 or (slen - send) > 10): #keeping a margin of 10
        difference1 = (0 - sstart) * 3 #this value will be negative
        difference2 = (slen - send) * 3
        new_start = start+difference1
        if new_start < 0:
            new_start = start
        new_stop = stop+difference2
        #print(start, stop, new_start, new_stop)
        extended_CDS = genome[new_start: new_stop]
        #print("new start, new stop, diff1 and diff 2:", new_start, new_stop, difference1, difference2, "\n")
        print("extended CDS\n", extended_CDS)
        print ("extended CDS length ", len(extended_CDS), "\n")
        new_df = pd.DataFrame(data = {'protein_name': name, 'CDS': extended_CDS, 'start': new_start, 'stop': new_stop}, index = [0])
        print ('new df after extension\n', new_df)
        return new_df

    #gotta address the time when sstart could be bigger than send, i.e. opposite matching
    elif sstart > send and ((send - 0) > 10 or (slen - sstart) > 10):
        difference2 = (send - 0) * 3 #this value should be positive
        difference1 = (sstart -slen) * 3 #this value shd be negative
        new_start = start+difference1
        if new_start < 0:
            new_start = start
        new_stop = stop+difference2
        #print("new start, new stop, diff1 and diff 2:", new_start, new_stop, difference1, difference2, "\n")
        extended_CDS = genome[new_start: new_stop]
        print("extended CDS\n", extended_CDS)
        print ("extended CDS length ", len(extended_CDS), "\n")
        new_df = pd.DataFrame(data = {'protein_name': name, 'CDS': extended_CDS, 'start': new_start, 'stop': new_stop}, index = [0])
        print ('new df after extension\n', new_df)
        return new_df


    else:
        return CDS_df.reset_index(drop=True)


#new blastx will be performed in the new extended CDS
