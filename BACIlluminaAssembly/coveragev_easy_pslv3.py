import os, sys, traceback

def proc_argv(argv=sys.argv):
  
    usage = """\

Note: Script for statistic work of the blat result. Specially used for correct contigs summary. 

Usage: *.py [-match num(default 0.98)] [-minlength num(default 0)] [-outmulti outmulti] [-outsingle outsingle] -query query.fa -target target.fa -psl psl_file -out output -depth depthfile"""

    match=0.98
    minlength=0
    outmulti=''
    outsingle=''
    
    if len(argv)>=9:
        num=len(argv)
        i = 1
        while i < len(argv) - 1:
            if argv[i] == "-match":
                i = i + 1
                match=float(argv[i])
                num=num-1
            elif argv[i] == "-minlength":
                i = i + 1
                minlength=float(argv[i])
                num=num-1
            elif argv[i] == "-psl":
                i = i + 1
                psl_file=(argv[i])
                num=num-1
            elif argv[i] == "-out":
                i = i + 1
                output=(argv[i])
                num=num-1
            elif argv[i] == "-outmulti":
                i = i + 1
                outmulti=(argv[i])
                num=num-1
            elif argv[i] == "-outsingle":
                i = i + 1
                outsingle=(argv[i])
                num=num-1
            elif argv[i] == "-query":
                i = i + 1
                query_file=(argv[i])
                num=num-1
            elif argv[i] == "-target":
                i = i + 1
                target_file=(argv[i])
                num=num-1
            elif argv[i] == "-depth":
                i = i + 1
                depthfile=(argv[i])
                num=num-1
            i = i + 1
    
    if len(argv) < 5 or num < 3:
        print(usage)
        sys.exit(1)
    assert match<=1,"invalidated input for match"

    
  
    assert os.path.exists(psl_file),"psl executable not found"
    return depthfile,match,minlength,psl_file,output,outmulti,outsingle,query_file,target_file
   

def pslfilter(psl,match,minlength,outmulti,outsingle):
    # print "psl filter..."
    title=[]
    match_multi={}
    match_uniq={}
    record_multi={}
    record_uniq={}
    name=[]
    outlinemulti=[]
    outlinesingle=[]
    count=0
    for eachline in psl:
        count+=1
        writeline=eachline
        try:
            eachline=(eachline.strip()).split()
            float(eachline[0])
        except:
            title.append(writeline)
            continue
            
            

        #if ((float(eachline[0])+float(eachline[3]))/float(eachline[10]))>float(match) and float(eachline[10])>minlength:  shut down filter
        if True:
            if eachline[17]=='1':
                if outmulti!='':
                    outlinemulti.append(writeline) 
                match_multi[eachline[9]]=float(eachline[10])
                blockSizes=int((eachline[18]).split(",")[0])
                try:
                    for i in range(int(eachline[15]),blockSizes+int(eachline[15])):
                        record_multi[eachline[13]][i]+=1                                         
                except:
                    #count+=1
                    record_multi[eachline[13]]=[0 for c in range(int(eachline[14]))]
                    for i in range(int(eachline[15]),blockSizes+int(eachline[15])):
                        record_multi[eachline[13]][i]+=1
                if eachline[9] not in name:
                    if outsingle!='':
                        outlinesingle.append(writeline)
                    match_uniq[eachline[9]]=float(eachline[0])
                    name.append(eachline[9])
                    try:
                        for i in range(int(eachline[15]),blockSizes+int(eachline[15])):
                            record_uniq[eachline[13]][i]+=1                                         
                    except:
                        #count+=1
                        record_uniq[eachline[13]]=[0 for c in range(int(eachline[14]))]
                        for i in range(int(eachline[15]),blockSizes+int(eachline[15])):
                            record_uniq[eachline[13]][i]+=1
                else:
                    try:
                        
                        outlinesingle.remove(writeline)
                        del match_uniq[eachline[9]]
                    except:
                        pass
            else:
                tends=[]
                qends=[]
                #tgaps=[]
                #qgaps[]
                gapMinus=[]
                canwrite=1
                tStarts=(eachline[20]).split(",")
                blockSizes=(eachline[18]).split(",")
                qStarts=(eachline[19]).split(",")
                blockCount=int(eachline[17])
                #print blockCount
                for i in range(blockCount):
                    #print i
                    qend=int(qStarts[i])+int(blockSizes[i])-1
                    tend=int(tStarts[i])+int(blockSizes[i])-1
                    try:
                        qgap=int(qStarts[i+1])-qend
                        tgap=int(tStarts[i+1])-tend
                    except:
                        pass
                    
                    tends.append(tend)
                    qends.append(qend)
                    #tgaps.append(tgap)
                    #qgaps.append(qgap)
                    gapMinus=abs(tgap-qgap)
                    if gapMinus>20:
                        canwrite=0
                        continue
                    try:
                        if (int(tStarts[i])<int(tEnds[i-1])) or (int(qStarts[i])<int(qends[i-1])):
                            canwrite=0
                            continue
                    except:
                        pass
                #print canwrite
                if canwrite==1 or canwrite==0:#shut down filter of gap
                    if outmulti!='':
                        outlinemulti.append(writeline) 
                    match_multi[eachline[9]]=float(eachline[0])
                    #long=len(tStarts)
                    for i in range(blockCount):
                        try:
                            for m in range(int(tStarts[i]),int(tends[i])):
                                record_multi[eachline[13]][m]+=1                                         
                        except:
                        #count+=1
                            record_multi[eachline[13]]=[0 for c in range(int(eachline[14]))]
                            #print tStarts,tends
                            for m in range(int(tStarts[i]),int(tends[i])):
                                record_multi[eachline[13]][m]+=1
                    if eachline[9] not in name:
                        if outsingle!='':
                            outlinesingle.append(writeline)
                        match_uniq[eachline[9]]=float(eachline[0])
                        name.append(eachline[9])
                        for i in range(blockCount):
                            try:
                                for m in range(int(tStarts[i]),int(tends[i])):
                                    record_uniq[eachline[13]][m]+=1
                            except:
                                record_uniq[eachline[13]]=[0 for c in range(int(eachline[14]))]
                                for m in range(int(tStarts[i]),int(tends[i])):
                                    record_uniq[eachline[13]][m]+=1
                            

                    else:
                        try:
                            del match_uniq[eachline[9]]
                            outlinesingle.remove(writeline)
                        except:
                            pass
                    
    # print "Done filtering the psl file\n"
    return title,match_multi,match_uniq,record_multi,record_uniq,outlinemulti
            

            
def faread(file):
    count=0
    length=0
    for eachline in file:
        if ">" in eachline :
            count+=1
        else:
            length+=len(eachline.strip())
    return count,length
        
def calcov(record):
    matchlong=0
    for eachkey in record:
        for position in record[eachkey]:
            if position>=1:
                matchlong+=1
    return matchlong

def caldict(match_dict):
    count=0
    length=0
    for key in match_dict:
        count+=1
        length+=float(match_dict[key])
    return count,length
    
 
def main():
    depth,match,minlength,psl_file,output_file,outmulti,outsingle,query_file,target_file = proc_argv()
    #print match,mis,minlength,psl_file,ref_file,output

    psl=open(psl_file,'r')
    query=open(query_file,'r')
    target=open(target_file,'r')
    
    query_count,query_length=faread(query)
    # print "query reading done..." 
    target_count,target_length=faread(target)
    # print "target reading done...\n" 
    query.close()
    target.close()

    title,match_multi,match_uniq,record_multi,record_uniq,outlinemulti=pslfilter(psl,match,minlength,outmulti,outsingle)

    match_count_multi,match_length_multi=caldict(match_multi)
    match_count_uniq,match_length_uniq=caldict(match_uniq)
    # print "cal coverage..."
    cov_muli=calcov(record_multi)
    cov_uniq=calcov(record_uniq)

    output=open(output_file,'w')
    # print 
    # output.write("Unique Matches Count:\t\t%d/%d(%lf)\n" % (match_count_uniq,query_count,(float(match_count_uniq)/float(query_count))))
    # output.write("Unique Matches Length:\t\t%d/%d(%lf)\n" % (match_length_uniq,query_length,(float(match_length_uniq)/float(query_length))))
    # output.write("All matches Uniquesites(Genome Coverage):\t\t%d/%d(%lf)\n\n" % (cov_uniq,target_length,(float(cov_uniq,)/float(target_length))))
    
    # output.write("Multiple Matches Count:\t\t%d/%d(%lf)\n" % (match_count_multi,query_count,(float(match_count_multi)/float(query_count))))
    # output.write("Multiple Matches Length:\t\t%d/%d(%lf)\n" % (match_length_multi,query_length,(float(match_length_multi)/float(query_length))))
    # output.write("All matches Multisites(Genome Coverage):\t\t%d/%d(%lf)\n\n" % (cov_muli,target_length,(float(cov_muli)/float(target_length))))

    output.close()
    psl.close()
    # print ("Unique Matches Count:\t\t%d/%d(%lf)\n" % (match_count_uniq,query_count,(float(match_count_uniq)/float(query_count))))
    # print ("Unique Matches Length:\t\t%d/%d(%lf)\n" % (match_length_uniq,query_length,(float(match_length_uniq)/float(query_length))))
    # print ("All matches Uniquesites(Genome Coverage):\t\t%d/%d(%lf)\n\n" % (cov_uniq,target_length,(float(cov_uniq,)/float(target_length))))
    
    # print ("Multiple Matches Count:\t\t%d/%d(%lf)\n" % (match_count_multi,query_count,(float(match_count_multi)/float(query_count))))
    # print ("Multiple Matches Length:\t\t%d/%d(%lf)\n" % (match_length_multi,query_length,(float(match_length_multi)/float(query_length))))
    # print ("All matches Multisites(Genome Coverage):\t\t%d/%d(%lf)\n\n" % (cov_muli,target_length,(float(cov_muli)/float(target_length))))
    
    
    #print readlong
    if outmulti!='':
        outm=open(outmulti,'w')
        for eachline in title:
            outm.write(eachline)
        for eachline in outlinemulti:
            outm.write(eachline)
        outm.close()

            
 
    psl.close()

    output.close()
    depthfile=open(depth,'w')
    for eachkey in record_multi:
        depthfile.write(">"+eachkey+'\n')
        count=0
        allcount=0
        for position in record_multi[eachkey]:
            depthfile.write(str(position)+' ')
            allcount+=1
            if position==0:
                count+=1
        depthfile.write('\n')
        # print ("    zero postion for %s: %d" % (eachkey,count))
        noncount=allcount-count
        # print ("non-zero postion for %s: %d" % (eachkey,noncount))
        # print ("     all postion for %s: %d" % (eachkey,allcount))
        percent=float(noncount)/float(allcount)
        print ("%s\t%d\t%d\t%.2f" % (eachkey,noncount,allcount,percent))
        # print ("%s\t%.2f" % (eachkey,percent))
        
            


main()