import sys 
import os

def get_precursor_features(Dir,whole_mature,whole_ID_mature,whole_ID, whole_miBase,\
                  whole_RF,whole_miBase_RF,fw_name):
    files=os.listdir(Dir) 
    fw=open("precursor_"+fw_name+"_tmp","w+")
    num_set=set()
    max_feature_ID_number=-1000
    min_feature_ID_numbe=1000
    max_feature_everyline_avg_in_block=-1000
    min_feature_everyline_avg_in_block=1000
    for file in files:
        Data_list=[]#save one file's data in a list
        reads_sum_in_file=0.0 # save the reads of one file 
        fr=open(Dir+"/"+file,"r")
        #first read one file , and the the sum of reads and process the data
        for line in fr.readlines()[1:]:
            s = line.strip().split("\t")
            reads_sum_in_file=reads_sum_in_file+float(s[2])
            Data_list.append(s) # the Data_list is the [[] ,[] .[]]
        fr.close()
        Data_list.sort(key=lambda x:(x[0],x[3]))#sort by mature and then ID 
        fr.close()
        #now start to iter the list
        pre_ID = Data_list[0][0] 
        dict_ID_rates={}  # save the mature-rates in dict
        dict_ID_reads={}  # save the mature_reads in dict
        dict_ID_mature_rates={} # save ID+mature-rates in dict
        dict_matures={}  # save mature-ID in dict
        block_list=[]  # save the same mature block 
        reads_sum_in_block=0.0
        feature_everyline_avg_in_block = 0.0
        #print "the Data_list len is "+str(len(Data_list))
        for item in Data_list:
            cur_ID = item[0]
            if cur_ID != pre_ID:
               #process the block list
               # the key is pre_mature 
               dict_ID_reads[pre_ID]=reads_sum_in_block
               dict_ID_rates[pre_ID]=reads_sum_in_block/reads_sum_in_file
               pre_mature = block_list[0][3]
               count_matures = 1 
               reads_sum_in_ID_mature_block=0.0
               for item2 in block_list:
                   #print item2
                   cur_mature = item2[3]
                   if pre_mature != cur_mature:
                      count_matures+=1
                      key = pre_ID+"_"+pre_mature
                      #print key , str(reads_sum_in_ID_mature_block) , str(reads_sum_in_block)
                      # save the different ID_mature\s rates in the dict
                      dict_ID_mature_rates[key]=reads_sum_in_ID_mature_block/reads_sum_in_block
                      pre_mature = cur_mature 
                      reads_sum_in_ID_mature_block = float(item2[2])
                   else:
                      reads_sum_in_ID_mature_block+=float(item2[2])
               key = cur_ID+"_"+pre_mature
               #print key , str(reads_sum_in_ID_mature_block) , str(reads_sum_in_block)
               dict_ID_mature_rates[key]=reads_sum_in_ID_mature_block/reads_sum_in_block
               dict_matures[pre_ID]=count_matures
               feature_everyline_avg_in_block+=1/float(len(block_list))
               pre_ID=cur_ID
               block_list=[]
               block_list.append(item)
               reads_sum_in_block = float(item[2])      
            else:
               reads_sum_in_block+=float(item[2])
               block_list.append(item) 
        dict_ID_reads[cur_ID]=reads_sum_in_block
        dict_ID_rates[cur_ID]=reads_sum_in_block/reads_sum_in_file
        pre_mature=block_list[0][3]
        count_matures=1
        reads_sum_in_ID_mature_block=0.0
        for item2 in block_list:
             cur_mature = item2[3]
             if pre_mature!= cur_mature:
                count_matures+=1
                key = pre_ID+"_"+cur_mature
                dict_ID_mature_rates[key]=reads_sum_in_ID_mature_block/reads_sum_in_block
                pre_mature=cur_mature
                reads_sum_in_ID_mature_block = float(item2[2])
             else:
                reads_sum_in_ID_mature_block = float(item2[2])
        key = cur_ID+"_"+cur_mature
        dict_ID_mature_rates[key] = reads_sum_in_ID_mature_block/reads_sum_in_block
        dict_matures[cur_ID]=count_matures
        feature_everyline_avg_in_block+=1/float(len(block_list))
        feature_everyline_avg_in_block/=float(len(dict_ID_reads))
        # end of the dict feature create 
        feature_ID_number=len(dict_ID_reads)
        feature_line=""
        feature_line=str(feature_ID_number)+"\t"+str(feature_everyline_avg_in_block)
        if feature_ID_number > max_feature_ID_number:  # get the max and min for normal
            max_feature_ID_number = feature_ID_number
        if feature_ID_number < min_feature_ID_numbe:
            min_feature_ID_numbe = feature_ID_number
        if feature_everyline_avg_in_block > max_feature_everyline_avg_in_block:
            max_feature_everyline_avg_in_block = feature_everyline_avg_in_block
        if feature_everyline_avg_in_block < min_feature_everyline_avg_in_block:
            min_feature_everyline_avg_in_block = feature_everyline_avg_in_block
        list_ID_rates=[]
        list_ID_matures=[]
        list_mature_ID_mature_rates=[]
        # get the feature list of mature-rates
        for item in whole_ID:
            if item in dict_ID_rates.keys():
                 tmp_list=[]
                 tmp_list.append(item)
                 tmp_list.append(dict_ID_rates[item])
                 list_ID_rates.append(tmp_list)
            else:
                 tmp_list=[]
                 tmp_list.append(item)
                 tmp_list.append(0)
                 list_ID_rates.append(tmp_list)
        list_ID_rates.sort(key=lambda x:(x[0]))
        for item in list_ID_rates:
            feature_line = feature_line+"\t"+str(item[1])
        #get the feature list of ID-mature rates 
        for item in whole_ID_mature:
            if item in dict_ID_mature_rates.keys():
                tmp_list=[]
                tmp_list.append(item)
                tmp_list.append(dict_ID_mature_rates[item])
                list_mature_ID_mature_rates.append(tmp_list)
            else:
                tmp_list=[]
                tmp_list.append(item)
                tmp_list.append(0)
                list_mature_ID_mature_rates.append(tmp_list)
        list_mature_ID_mature_rates.sort(key=lambda x:(x[0]))
        for item in list_mature_ID_mature_rates:
             feature_line = feature_line+"\t"+str(item[1])
        #get the feature list of ID-mature-rates
        for item in whole_ID:
            if item in dict_matures.keys():
                 tmp_list=[]
                 tmp_list.append(item)
                 tmp_list.append(dict_matures[item])
                 list_ID_matures.append(tmp_list)
            else:
                 tmp_list=[]
                 tmp_list.append(item)
                 tmp_list.append(0)
                 list_ID_matures.append(tmp_list)
        list_ID_matures.sort(key=lambda x:(x[0]))
        for item in list_ID_matures:
            if item[1]==0:
                feature_line = feature_line+"\t"+str(0)
            else:
                feature_line = feature_line+"\t"+str(1/float(item[1]))
        feature_line+="\n"
        fw.write(feature_line)
    fw.close()
    fr = open("precursor_"+fw_name+"_tmp","r")
    fw = open("precursor_"+fw_name,"w+")
    for line in fr: # nromal the feature 1 and feature 2
         s= line.strip().split("\t")
         if max_feature_ID_number!=min_feature_ID_numbe:
             s[0]=(float(s[0])-min_feature_ID_numbe)/(max_feature_ID_number-min_feature_ID_numbe)
         else:
             s[0]=1.0
         s[0]=str(s[0])
         if max_feature_everyline_avg_in_block!=min_feature_everyline_avg_in_block:
             s[1]=(float(s[1])-min_feature_everyline_avg_in_block)/(max_feature_everyline_avg_in_block-min_feature_everyline_avg_in_block)
         else:
             s[1]=1.0
         s[1]=str(s[1])
         feature_line="\t".join(s)+"\n"
         fw.write(feature_line)
    fr.close()
    fw.close()
    os.remove("precursor_"+fw_name+"_tmp")

#        print "the precursor feature number="+str(len(whole_ID))+"+"+str(len(whole_ID_mature))+"+"+\
#               str(len(whole_ID))+"+4 ="+str(len(whole_ID)+len(whole_ID_mature)+len(whole_ID)+2)
#        print "the real featuer number is : "+str(len(feature_line.strip().split("\t")))
#        if len(whole_ID)+len(whole_ID_mature)+len(whole_ID)+2!=len(feature_line.strip().split("\t")):
#            print "the precursor is no equal"
#            exit()
    #end 
 
