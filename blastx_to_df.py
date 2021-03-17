#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import re
import pandas as pd


# In[2]:


#in this function see if it is a single hit
def remove_duplicates(blastx_df):
    #first remove any duplicate matches, if qend of another match is smaller, then its a duplicate
    x = 0
    for i in range(0,len(blastx_df)):
        blastx_df = blastx_df.reset_index(drop = True)
        if i ==0:
            continue
        qend1 = blastx_df.loc[i-1-x, 'qend']
        #print ("qend 1 and 2 ", qend1)
        qend2 = blastx_df.loc[i-x,'qend']
        #print (qend2)
        if qend2 <= qend1:
            print ('removed', str(blastx_df.iloc[i-x,:]))
            blastx_df = blastx_df.drop([i-x])
            #the one we will check will have to stay same
            x = x+1
    blastx_df = blastx_df.reset_index(drop = True)
    print ('\n after removing duplicates\n', str(blastx_df))
    return blastx_df


def check_single_hit(blastx_df):
    length = len(blastx_df)
    if length == 1:
        #print ("single hit\n")
        return "yes"
    else:
        return "no"


#make a function to run blastx and make a result dataframe for the query CDS
#the return should be a dataframe with qstart, qend and all other values
def blastx_to_df(CDS,db):
    print ("blastx db is ", db)
    CDS1 = open('CDS', 'w')
    CDS1.write('>query_CDS\n')
    CDS1.write(CDS)
    CDS1.close()
    cmd = "blastx -query CDS -db %s -max_target_seqs 1 -outfmt '6 qseqid sseqid qstart qend sstart send qlen slen qframe sframe qseq sseq stitle pident' > result_trial" %(db)
    os.system(cmd)
    if os.path.getsize("result_trial") > 0:
        result = pd.read_csv('result_trial', sep= "\t", header = None)
        #print ("result: ", result)
    else:
        return 0 #if no hits, 0 is returned
    length = len(result)
    num_hits = 0
    sid_prev = 0
    for i in range(length): #got the number of hits for a CDS
        sseqid = result.iloc[i, 1]
        if sid_prev == sseqid:
            #text = 'it has more than 1 hits to sequence id \n'+ str(sseqid)+ '\n'
            #print(text)
            num_hits += 1
        elif num_hits > 0 and sid_prev != sseqid:
            break
        else:
            num_hits += 1
            sid_prev = sseqid

    columns = ['qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'stitle', 'pident', 'qseq']
    print('num of hits %s' %str(num_hits))

    #declare all the lists
    qstart = [0]* num_hits
    qend = [0]* num_hits
    sstart = [0]* num_hits
    send = [0]* num_hits
    qframe = [0] * num_hits
    qlen = [0] * num_hits
    slen = [0] * num_hits
    stitle = [0] * num_hits
    pident = [0] * num_hits
    qseq = [0] * num_hits
    stats= []

    #get the qstart and qend positions
    for i in range(num_hits):
        qstart[i] = result.iloc[i, 2] -1 #-1 to account for py index
        qend[i] = result.iloc[i,3] - 1
        sstart[i] = result.iloc[i, 4] -1
        send[i] = result.iloc[i,5]-1
        qframe[i] = result.iloc[i,8]
        qlen[i] = result.iloc[i,6]
        slen[i] = result.iloc[i,7]
        stitle[i] = result.iloc[i,12]
        pident[i] = result.iloc[i,13]
        qseq[i] = result.iloc[i,10]
        stats.append([qstart[i], qend[i], sstart[i], send[i], qlen[i], slen[i], stitle[i], pident[i], qseq[i]])
        #print (stats)

    #avoid the wrongly concatenated sequences from being corrected
    #for i in range(1,num_hits): #so that repeated CDSes are not confused with split ones
        #skips if the CDSes give same frame reading
        #if qframe[i] == qframe[i-1]:
            #report.write("But its a repeat not split protein sequence since the qframe are same for two matches \n")
            #return CDS

    #fix if reverse matching
    if qend[0] < qstart[0]:
        stats = []
        for i in range(num_hits):
            temp = qstart[i]
            qstart[i] = qend[i]
            qend[i] = temp
            temp2 = sstart[i]
            sstart[i] = send[i]
            send[i] = temp2
            stats.append([qstart[i], qend[i], sstart[i], send[i], qlen[i], slen[i], stitle[i], pident[i], qseq[i]])

    #build a dataframe for the qstart and qend positions
    df = pd.DataFrame(data= stats, columns = columns)
    df = df.sort_values(by = "qstart", ignore_index=True)
    print("the blastx dataframe: \n", str(df))
    #print ('return from blastx_to_df function\n', df)
    return df


# In[3]:

#return only when single hit
