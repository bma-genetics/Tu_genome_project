#!/usr/bin/env python
import sys
import collections

def read_pair():
    pair_dic = {}
    infile = open('/mnt/Data_diskD/Wheat_Agenome_DNA/repeat_TE_chr/LTR_homolog_search/paired_LTR_search.out', 'r')   
    for line in infile.readlines():
        if not line.strip():continue 
        try:
            line = line.rstrip()
            word = line.split()
            pair_dic[word[9][1:]] = ""
            print word[9][1:]+"**"
        except ValueError:
            continue
    infile.close()
    return pair_dic

def read_TE(pair_dic1):
    TE_dic = {}
    infile = open('/mnt/Data_diskD/Wheat_Agenome_DNA/repeat_TE_chr/LTR_homolog_search/Tu_classI_TE.fasta', 'r')   
    for line in infile.readlines():
        if not line.strip():continue 
        try:
            if(">" in line):
                line = line.rstrip()
                word = line.split()
                if pair_dic1.has_key(word[1]): 
                    del pair_dic1[word[1]]
                    print word[0][4:]+"***"+word[1]+"***"+word[2][10:]+"***"+word[3]+"***"+word[4]
                    TE_dic[int(word[0][4:])] = TE_dic.get(int(word[0][4:]),{})
                    TE_dic[int(word[0][4:])][int(word[3])] = TE_dic[int(word[0][4:])].get(int(word[3]),{})
                    TE_dic[int(word[0][4:])][int(word[3])][int(word[4])]=TE_dic[int(word[0][4:])][int(word[3])].get(int(word[4]),word[2][10:])
                    #TE_dic[int(word[0][4:])][int(word[3])][int(word[4])]=word[2][10:]
        except ValueError:
            continue
    infile.close()
    return TE_dic


def iden_TE(TE_dic):
    #print len(TE_dic.keys())
    outfile = open("iden_gapsy_copia.out", "w")
    outfile.write("chr s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr type\n")
    for i in range(1,8,1):
        od = collections.OrderedDict(sorted(TE_dic[i].items()))
        print i
        infile = open('/mnt/Data_diskD/Wheat_Agenome_DNA/repeat_TE_chr/LTR_digest/ltrharvest_te_chr'+str(i), 'r')   
        for line in infile.readlines():
            if not line.strip():continue 
            try:
                if "#" in line:continue
                #print line
                line = line.rstrip()
                word = line.split()
                print word[0]+" "+word[1]+" "+word[2]
                for k, v in od.items():
                    od1 = collections.OrderedDict(sorted(v.items()))
                    for ik, iv in od1.items():
                        if(ik < int(word[0])): 
                            del od[k][ik]
                            continue
                        if(k > int(word[1])): break
                        x = range(int(word[0]),int(word[1]))
                        y = range(k,ik)
                        olap = range(max(x[0], y[0]), min(x[-1], y[-1])+1)
                        #print x
                        #print y
                        #print olap
                        if(len(olap) > 0.5*int(word[2])):
                           outfile.write(str(i)+"\t"+line+"\t")
                           outfile.write(iv)
                           outfile.write("\n")
            except ValueError:
                continue
        infile.close()
    outfile.close()

if __name__ == '__main__':
    pair_dic1 = read_pair()
    TE_dic = read_TE(pair_dic1)
    iden_TE(TE_dic)
    
    
