#iterative blastx for the changes in the overlapping nts
import re
import os
import sys
import pandas as pd
from blastx_to_df import remove_duplicates

#a program to fix the multiple frameshift of a CDS based on its hits against protein db
#currently limited to genes giving only two hits against a protein database
#could extend it to iteratively fix hits in pairs
#this has a problem: sometimes same CDS shows hits against different proteins when read at different frames



def blastx_cmd(CDS_containing_file, db):
    
    if CDS_containing_file == 'CDSfile':
        cmd = "blastx -query CDSfile -db %s -max_target_seqs 1 -soft_masking True -outfmt '6 qseqid sseqid \
        qstart qend sstart send qlen slen qframe sframe qseq sseq stitle pident bitscore' > result_trial" %(db)
    
    else:
        new_CDS_file = open('new_CDS_file', 'w+')
        new_CDS_file.write('>new_CDS\n' + CDS_containing_file)
        new_CDS_file.close()
        cmd = "blastx -query new_CDS_file -db %s -max_target_seqs 1 -soft_masking True -outfmt '6 qseqid sseqid \
        qstart qend sstart send qlen slen qframe sframe qseq sseq stitle pident bitscore' > result_trial" %(db)

    print ('cmd is', cmd)
    
    os.system(cmd)
    if os.path.getsize("result_trial") > 0:
        result = pd.read_csv('result_trial', sep= "\t", header = None)
    else:
        print ('no hits to protein db')
        return 0

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

    columns = ['qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'stitle', 'qseq', 'bitscore']
    #print('num of hits %s' %str(num_hits))

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
    bitscore = [0] * num_hits
    stats= []


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
        bitscore[i] = result.iloc[i,14 ]
        stats.append([qstart[i], qend[i], sstart[i], send[i], qlen[i], slen[i], stitle[i], bitscore[i], qseq[i]])

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
            stats.append([qstart[i], qend[i], sstart[i], send[i], qlen[i], slen[i], stitle[i], bitscore[i], qseq[i]])

    #build a dataframe for the qstart and qend positions
    df = pd.DataFrame(data= stats, columns = columns)
    df = df.sort_values(by = "qstart", ignore_index=True)
    print("the blastx dataframe: \n", df)


    blastx_df = remove_duplicates(df)
    num_hits = len(blastx_df)
    df_and_hits = [blastx_df, num_hits]
    return df_and_hits

#-----------------------------------------------------------------------------------------------------------------------------


# def remove_duplicates(blastx_df):
#     #first remove any duplicate matches, if qend of another match is smaller, then its a duplicate
#     x = 0
#     for i in range(0,len(blastx_df)):
#         blastx_df = blastx_df.reset_index(drop = True)
#         if i ==0:
#             continue
#         qend1 = blastx_df.loc[i-1-x, 'qend']
#         #print ("qend 1 and 2 ", qend1)
#         qend2 = blastx_df.loc[i-x,'qend']
#         #print (qend2)
#         if qend2 <= qend1:
#             print ('removed', blastx_df.iloc[i-x,:])
#             blastx_df = blastx_df.drop([i-x])
#             #the one we will check will have to stay same
#             x = x+1
#     blastx_df = blastx_df.reset_index(drop = True)
#     print ('\n after removing duplicates\n', blastx_df)
#     return blastx_df


#------------------------------------------------------------------------------------------------------------------------
def find_bitscore_blastn(CDS, nt_db):
    new_CDS_file = open("CDS2", 'w+')
    new_CDS_file.write('>new_CDS\n')
    new_CDS_file.write(CDS+'\n')
    new_CDS_file.close()
    cmd = "blastn -query CDS2 -db %s -max_target_seqs 1 -perc_identity 80 -soft_masking True -outfmt '6 bitscore' > res" %(nt_db)
    os.system(cmd)
    #print('cmd is', cmd)
    bitscore = 0
    if os.path.getsize("res") > 0:
        result = pd.read_csv('res', sep= "\t", header = None)
        bitscore = result.iloc[0,0]
    return bitscore

def find_bitscore_blastx(CDS, db):
    new_CDS_file = open("CDS2", 'w+')
    new_CDS_file.write('>new_CDS\n')
    new_CDS_file.write(CDS+'\n')
    new_CDS_file.close()
    cmd = "blastx -query CDS2 -db %s -max_target_seqs 1 -soft_masking True -outfmt '6 bitscore' > res" %(db)
    os.system(cmd)
    #print('cmd is', cmd)
    bitscore = 0
    if os.path.getsize("res") > 0:
        result = pd.read_csv('res', sep= "\t", header = None)
        bitscore = result.iloc[0,0]
    return bitscore

#---------------------------------------------------------------------------------------------------------------------------
def fix_insertion(CDS, endpoint, check_region_len, db, nt_db):
    print ('check_region_len', check_region_len)
    new_CDS_file = open("CDS2", 'w+')
    new_CDS_file.write('>new_CDS\n')
    new_CDS_file.write(CDS+'\n')
    new_CDS_file.close()
    blastn_cmd = "blastn -query CDS2 -db %s -max_target_seqs 1 -perc_identity 80 -soft_masking True -outfmt '6 bitscore' > res" %(nt_db)
    os.system(blastn_cmd)
    blastn = False
    #whether use blastn or blastx for changing the sequence
    if os.path.getsize("res") > 0:
        print ('using blastn for fixing insertion')
        blastn = True
    blastn = False #cause in this case I don't want to use blastn again for fixing. 
    print("check region: ", CDS[endpoint-24:endpoint+24])
    bitscores = []
    js = []
    new_CDSes = []
    for j in range(check_region_len): #need to throw out a nt
        CDS_parts_change = [0] * 2
        #print (' It has insertion and j is ', j)
        # print (CDS[endpoint+j-2:endpoint+j+2])
        # print(CDS[endpoint+j-1:endpoint+j+3])
        if not (re.search(r'[A]{2}[CGT]|[C]{2}[ATG]|[G]{2}[ACT]|[T]{2}[ACG]', CDS[endpoint+j-1:endpoint+j+2]) or re.search(r'[CGT][A]{2}|[ATG][C]{2}|[ACT][G]{2}|[ACG][T]{2}', CDS[endpoint+j-1:endpoint+j+2])):
            continue

        #it will remove the j only if its one of the homopolymers

        print("homopolymeric region: ", CDS[endpoint+j-2:endpoint+j+3])
        CDS_parts_change[0] = CDS[:endpoint] + CDS[endpoint:endpoint+j]  + CDS[endpoint+j+1:endpoint+check_region_len] #iteratively remove each 'nt'
        CDS_parts_change[0+1] = CDS[endpoint + check_region_len:]
        # print(CDS_parts_change[i])
        # print(CDS_parts_change[i+1])
        new_CDS = CDS_parts_change[0] + CDS_parts_change[0+1]
        print ('j - value',j, '\n')
        if blastn == True:
            bitscore = find_bitscore_blastn(new_CDS, nt_db)
        else:
            bitscore = find_bitscore_blastx(new_CDS, db)
        print ('bitscore', bitscore, '\n\n')
        if bitscore>0:
            bitscores.append(bitscore)
            js.append(j)
            new_CDSes.append(new_CDS)
    if len(bitscores) > 0:
        max_bit_index = bitscores.index(max(bitscores))
        j_value = js[max_bit_index]
        #print ('max j_value ',j_value )
        new_CDS = new_CDSes[max_bit_index]
        #print (new_CDS)
        change_pos_CDS = endpoint + j_value 
    else:
        change_pos_CDS = 0
        new_CDS = CDS
    list1 = [change_pos_CDS, new_CDS]
    return list1



#-----------------------------------------------------------------------------------------------------------------------------


def fix_deletion(CDS, endpoint, check_region_len, db, nt_db):
    new_CDS_file = open("CDS2", 'w+')
    new_CDS_file.write('>new_CDS\n')
    new_CDS_file.write(CDS+'\n')
    new_CDS_file.close()
    blastn_cmd = "blastn -query CDS2 -db %s -max_target_seqs 1 -perc_identity 80 -soft_masking True -outfmt '6 bitscore' > res" %(nt_db)
    os.system(blastn_cmd)
    blastn = False
    os.path.getsize("res") 
    #whether use blastn or blastx for changing the sequence
    if os.path.getsize("res") > 0:
         print ("using blastn for fixing deletion")
         blastn = True
    blastn = False
    print ('check_region_len', check_region_len)
    print("check region: ", CDS[endpoint-24:endpoint+24])
    bitscores = []
    new_CDSes = []
    js = []
    change_pos_CDS = 0
    new_CDS = CDS
    for j in range(check_region_len):
        #print("check region: ", CDS[endpoint+j-2:endpoint+j+3])
        #print (' It has deletion and j is ', j)
        CDS_parts_change = [0] * 2
        print ('j ',j )
        #print (CDS[endpoint+j-2:endpoint+j+1])
        #print(CDS[endpoint+j-1:endpoint+j+2])
        if not (re.search(r'[A]{1}[CGT]|[C]{1}[AGT]|[G]{1}[CAT]|[T]{1}[CGA]', CDS[endpoint+j-1:endpoint+j+1])):
             continue
        else:
            fix_with = CDS[endpoint+j]
        print("homopolymeric region: ", CDS[endpoint+j-2:endpoint+j+3])
        CDS_parts_change[0] = CDS[:endpoint] + CDS[endpoint:endpoint+j]+ fix_with + CDS[endpoint+j:endpoint+check_region_len] #iteratively insert  each 'nt'
        CDS_parts_change[0+1] = CDS[endpoint + check_region_len:]
        new_CDS = CDS_parts_change[0] + CDS_parts_change[0+1]
        #print ('j - vale and cds ',j, '\n', new_CDS)
        if blastn == True:
            bitscore = find_bitscore_blastn(new_CDS, nt_db)
        else:
            bitscore = find_bitscore_blastx(new_CDS, db)
        print ('bitscore', bitscore, '\n\n')
        bitscores.append(bitscore)
        js.append(j)
        new_CDSes.append(new_CDS)
    if len(bitscores) > 0:
        max_bit_index = bitscores.index(max(bitscores))
        j_value = js[max_bit_index]
        #print ('max j_value ',j_value )
        new_CDS = new_CDSes[max_bit_index]
        #print (new_CDS, '\n')
        change_pos_CDS = endpoint + j_value 
        #print (change_pos_CDS)
    else:
        change_pos_CDS = 0
        new_CDS = CDS
        #print (change_pos_CDS)
    list1 = [change_pos_CDS, new_CDS]
    return list1


#-----------------------------------------------------------------------------------------------------------------------------    

def blastx_CDS(CDSfile, CDSes_fixed, blastx_df, CDS_start, nt_shift, db, nt_db):
    new_CDS = 0
    
    CDS_read = open(CDSfile, 'r')
    CDS = ''
    # blastxDF_and_hits = blastx_cmd(CDSfile,db)
    # if isinstance(blastxDF_and_hits, int):
    #     return 0
    for i in CDS_read.readlines():
        if not re.search('>', i):
            CDS = CDS + i
        else:
            header = i
    # blastx_df = blastxDF_and_hits[0]
    # num_hits = blastxDF_and_hits[1]
    #CDS.replace('-', '')
    num_hits = len(blastx_df)
    #print ('CDS is', CDS)
    #CDS_parts = [0] * num_hits
    #CDS_parts_change = [0] * num_hits
    insertion = False
    deletion = False
    qstart = [0]* num_hits
    qend = [0]* num_hits
    sstart = [0]* num_hits
    send = [0]* num_hits
    for i in range(num_hits):
        qstart[i] = blastx_df.loc[i,'qstart']
        qend[i] = blastx_df.loc[i,'qend']
        sstart[i] = blastx_df.loc[i,'sstart']
        send[i] = blastx_df.loc[i,'send']

    for i in range(num_hits-1):
        print (i, qend[i], qstart[i])
        if qend[i] >= qstart[i+1]:
            overlap = qend[i] - qstart[i+1] + 1
            check_region_len = overlap + 48 #24 nt up and downstream of the overlap sites?
            print ('overlap', overlap)
            
            endpoint = qstart[i+1] -24
            if overlap > 200:
                print ('overlap exceeds 200 nts. It should not be a frameshift error. something else is going on')
                CDSes_fixed.write(header+'overlap exceeds 200 nts. It should not be a frameshift error. something else is going on')
                return 0
                
            #print ('check_region len ', check_region_len, 'from', endpoint, 'to', qend[i] + 6 )
            #CDS_parts[i] = CDS[: endpoint]
            #x = len(CDS[qstart[0]:endpoint])
            #insertion or deletion?
            if (len(CDS[qstart[0]:endpoint]) ) %3 == 1:
                insertion = True

            elif (len(CDS[qstart[0]:endpoint])) % 3 == 2:
                deletion = True

            if insertion:
                print ('insertion found in the CDS. Will throw out a nt now')
                change_pos_CDS, new_CDS = fix_insertion(CDS, endpoint, check_region_len, db, nt_db)
                text = 'deleted %s in the string %s at 9th from start resulting in %s' %(CDS[change_pos_CDS], CDS[change_pos_CDS-8 : change_pos_CDS + 8], new_CDS[change_pos_CDS-8 : change_pos_CDS + 8])
                print (text)
                change_dict = {'change type': 'deletion', 'position in new  genome': CDS_start+nt_shift+change_pos_CDS, 'before': CDS[change_pos_CDS], 'after': '', 'seq before': CDS[change_pos_CDS-5:change_pos_CDS+6], 'seq after': new_CDS[change_pos_CDS-5:change_pos_CDS+6] }
                print ('change_dict')
            
            if deletion:
                print ('deletion found in the CDS. Will include a nt now')
                change_pos_CDS, new_CDS = fix_deletion(CDS, endpoint, check_region_len, db, nt_db)
                text = 'inserted %s in the string %s resulting in %s' %(new_CDS[change_pos_CDS], CDS[change_pos_CDS-8 : change_pos_CDS + 8], new_CDS[change_pos_CDS-8 : change_pos_CDS + 8])
                print (text)
                change_dict = {'change type': 'insertion', 'position in new genome': CDS_start+nt_shift+change_pos_CDS, 'before': '', 'after': new_CDS[change_pos_CDS], 'seq before': CDS[change_pos_CDS-5:change_pos_CDS+6], 'seq after': new_CDS[change_pos_CDS-5:change_pos_CDS+6] }
                print ('change_dict')
 

        elif qend[i] < qstart[i+1]:
           
            missed_nts = qstart[i+1] - qend[i] - 1
            print("missed_nts ", missed_nts)
            #print ("it has missed nucleotides not translated to any protein")
            check_region_len = missed_nts + 48 #3 nts downstream and upstream of the first match ending site
            endpoint = qend[i] -24
            #print ('check_region len ', check_region_len, 'from', endpoint, 'to', qstart[1] + 6 )
          
            #CDS_parts[i] = CDS[: endpoint]
            if missed_nts %3 == 1:
                insertion = True
            elif missed_nts %3 ==2:
                deletion = True
            if insertion:
                print ('insertion found in the CDS. Will throw out a nt now')
                change_pos_CDS, new_CDS = fix_insertion(CDS, endpoint, check_region_len, db, nt_db)
                text = 'deleted %s in the string %s at 9th from start resulting in %s' %(CDS[change_pos_CDS], CDS[change_pos_CDS-8 : change_pos_CDS + 8], new_CDS[change_pos_CDS-8 : change_pos_CDS + 8])
                print (text)
                change_dict = {'change type': 'deletion', 'position in genome': CDS_start+nt_shift+change_pos_CDS, 'before': CDS[change_pos_CDS], 'after': '', 'seq before': CDS[change_pos_CDS-5:change_pos_CDS+6], 'seq after': new_CDS[change_pos_CDS-5:change_pos_CDS+6] }
                print ('change_dict')
            if deletion:
                print ('deletion found in the CDS. Will include a nt now')
                change_pos_CDS, new_CDS = fix_deletion(CDS, endpoint, check_region_len, db, nt_db)
                #print ('the fixed CDS is', new_CDS)
                text = 'inserted %s in the string %s resulting in %s' %(new_CDS[change_pos_CDS], CDS[change_pos_CDS-8 : change_pos_CDS + 8], new_CDS[change_pos_CDS-8 : change_pos_CDS + 8]) 
                print ('text')
                change_dict = {'change type': 'insertion', 'position in genome': CDS_start+nt_shift+change_pos_CDS, 'before': '', 'after': new_CDS[change_pos_CDS], 'seq before': CDS[change_pos_CDS-5:change_pos_CDS+6], 'seq after': new_CDS[change_pos_CDS-5:change_pos_CDS+6] }
                print ('change_dict')

    #last blastx to see if it really worked and produced a single hit
    
    if not new_CDS or not change_pos_CDS:
        CDSes_fixed.write(header + "Could not be fixed. No insertion or deletion found in the CDS \n")
        return 0
    
    changed_list = [new_CDS, change_dict]
    blastxDF_and_hits = blastx_cmd(new_CDS,db)
    if isinstance(blastxDF_and_hits, int):
        CDSes_fixed.write(header + "Could not be fixed. No hits in the protein db found \n")
    blastx_df = blastxDF_and_hits[0]
    num_hits = blastxDF_and_hits[1]
    
    if num_hits > 1:
        CDSes_fixed.write(header + "Could not be fixed. Due to multiple hits even before manual alteration \n")
        print ("Still gave multiple protein hits. Does not seem to be the result of an indel")
        return 0
    else:
        CDSes_fixed.write(header + new_CDS +  '\n' + text + '\n')
    return changed_list


#---------------------------------------------------------------------------------------------------------------------------

# def main():
#     #read a file containing CDSes
    
#     CDSes_file = open('CDSes_to_fix', 'r')
    
#     CDSes_fixed = open('CDSes_fixed_file', 'w+')
#     CDSes = []
#     CDS = ''
#     last_line = CDSes_file.readlines()[-1]
#     CDSes_file.seek(0)
#     for i in CDSes_file.readlines():
#         if i == last_line:
#             CDS = CDS.replace('-', '')
#             CDSes.append(CDS)
#         elif CDS and re.search('>', i):
#             CDS = CDS.replace('-', '')
#             CDSes.append(CDS)
#             CDS = ''
#             CDS = CDS+i
#         elif not CDS and re.search('>', i):
#             CDS = CDS+i
#             continue
#         else:
#             CDS = CDS + i.strip()
    
#     CDSes_file.close()
#     for each_CDS in CDSes:
#         CDSfile = open('CDSfile', 'w+')
#         CDSfile.write(each_CDS)
#         CDSfile.close()
#         blastx_CDS('CDSfile', CDSes_fixed)
#     CDSes_fixed.close()



# #-----------------------------------------------------------------------------------------------------------------------------

# if __name__ == "__main__":
#     main()

# #input: 
