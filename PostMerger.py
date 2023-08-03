import os, argparse
  
HLA=["HLA_A","HLA_B","HLA_C","HLA_DRB1","HLA_DQA1","HLA_DQB1","HLA_DPA1","HLA_DPB1"]
vcf_files=[]

def HLAextract_SNP2HLA(vcf_input):
    for i in HLA:
        os.system("grep '#CHROM' %s | sed -e 's/^#//' > %s.%s_tmp"%(vcf_input,vcf_input,i))
        os.system("grep '%s' %s >> %s.%s_tmp"%(i,vcf_input,vcf_input,i))

def VcfWeight(__inputlist__):
    vcf_name=__inputlist__[0]
    weight=float(__inputlist__[1])
    vcf=[]
    with open(vcf_name,"r") as f1:
        for i in f1: vcf.append(i.split())

    row=len(vcf)
    col=len(vcf[0])
    new_HLA_vcf=[[0 for i in range(col)] for j in range(row)]

    for i in range(row): new_HLA_vcf[i][:9]=vcf[i][:9]
    new_HLA_vcf[0][9:]=vcf[0][9:]
    
    for i in range(1,row):
        for j in range(9,col):
            HLApb=float(vcf[i][j].split(":")[2].split(",")[0])+float(vcf[i][j].split(":")[2].split(",")[1])/2
            new_HLA_vcf[i][j]=str(round(HLApb*weight,4))
            
    with open(vcf_name+"."+str(weight),"w") as f1:
        for i in new_HLA_vcf:
            f1.write("\t".join(i)+"\n")

def VcfMerge(HLA_allele):
    alleles_dict={}
    alleles_pp_dict={}
    alleles_merged=[]
    for i in range(len(vcf_files)):
        alleles_dict[i]={}
        with open(vcf_files[i][0]+".HLA_"+HLA_allele+"."+vcf_files[i][1],"r") as vcf:
            FID=next(vcf).split()[9:]
            for j in FID: alleles_dict[i][j]={}
            for j in vcf:
                pp=j.split()[9:]
                HLA_vcf=j.split()[2]
                for l in range(len(pp)):
                    alleles_dict[i][FID[l]][HLA_vcf]=pp[l]
        alleles_merged=alleles_merged+list(alleles_dict[i][FID[0]].keys())
    alleles_merged=list(set(alleles_merged))
    if not len(alleles_merged):
        with open(output+"."+HLA_allele+".alleles","w") as fw:
            fw.write("")
        return

    ### alleles_dict[vcf_files][FID][HLA] ###

    alleles_answer=[[0 for i in range(7)] for j in FID]
    cnt=0
    for i in FID:
        alleles_pp_dict[i]={}
        for j in alleles_merged:
            o,pp=0,0
            for k in range(len(vcf_files)):
                 if j in alleles_dict[k][i].keys():
                        o=o+float(vcf_files[k][1])
                        pp=pp+float(alleles_dict[k][i][j])
            if o!=0:
                alleles_pp_dict[i][j]=round(pp/o,4)
        alleles_answer[cnt][0:3]=i,i,HLA_allele
        alleles_answer[cnt][5]=max(alleles_pp_dict[i].values())
        allele_1=[k for k,v in alleles_pp_dict[i].items() if v==alleles_answer[cnt][5]][0]
        mmax,ta=0,True
        for k,v in alleles_pp_dict[i].items():
            if k!=allele_1 and v>=alleles_answer[cnt][5]/2 and mmax < v:
                mmax=v
                allele_2=k
                alleles_answer[cnt][6]=v
                ta=False
        if ta:
            allele_2=allele_1
            alleles_answer[cnt][6]=alleles_answer[cnt][5]

        alleles_answer[cnt][3]=allele_1.split("_")[-1][0:2]+","+allele_2.split("_")[-1][0:2]
        pp,qq=allele_1.split("_")[-1],allele_2.split("_")[-1]
        if HLA_allele=='DRB1':
            if allele_1.split("_")[-1]=='1454' : pp='1401'
            if allele_2.split("_")[-1]=='1454' : qq='1401'

        alleles_answer[cnt][4]=pp+","+qq
        cnt=cnt+1

    with open(output+"."+HLA_allele+".alleles","w") as fw:
        for i in alleles_answer:
            fw.write(" ".join(map(str,i))+"\n")

if __name__=="__main__":
    
    ### SNP2HLA ###
    parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("-i","--input", help = "\n VCF file list \n\n")
    parser.add_argument("-o","--output", help = "\n Output file prefix \n\n")
    args = parser.parse_args()

    vcf_list,output=args.input, args.output

    with open(vcf_list,"r") as f1:
        for i in f1:
            vcf_files.append(i.split())

    for i in vcf_files:
        HLAextract_SNP2HLA(i[0])
    
    for i in HLA:
        for j in vcf_files:
            with open(j[0]+"."+i+"_tmp","r") as f1 , open(j[0]+"."+i,"w") as f1_w:
                    f1_w.write(next(f1))
                    for k in f1:
                            if len(k.split()[2].split("_")[2]) >= 4 : f1_w.write(k)
            os.system("rm %s"%(j[0]+"."+i+"_tmp"))
    
    for i in HLA:
        for j in vcf_files:
            VcfWeight([j[0]+"."+i,j[1]])
        VcfMerge(i.split("HLA_")[1])
   
    alleles_list=""
    for i in HLA:
        alleles_list=alleles_list+output+"."+i.split("HLA_")[1]+".alleles"+"\t"
    
    os.system("cat %s > %s"%(alleles_list,output+".alleles"))
    for i in HLA: os.system("rm %s"%(output+"."+i.split("HLA_")[1]+".alleles"+"\t"))
    
    for i in HLA: 
        for j in vcf_files:
            os.system("rm %s"%(j[0]+"."+i+"*"))
    

