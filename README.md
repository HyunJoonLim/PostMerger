# PostMerger

## (1) Introduction

***PostMerger*** is a software package that merges multiple reference panels for HLA imputation. The merged reference panel could be used for running **SNP2HLA**.


*SNP2HLA* : <http://software.broadinstitute.org/mpg/snp2hla/>


---


## (2) Steps


***PostMerger*** consists of the following four steps.


The three steps of ***PostMerger*** for merging multiple imputation results into one.

**Step1**. Calculating the posterior probability.  
**Step2**. Combining poster probabilities at the given weights.  
**Step3**. Calling the best HLA allele markers.


***


## (3) Usage


### **Dependency**


* Python==3.7

Install **Python** (<http://www.python.org/downloads/>)   



### **RUN**

    USAGE: python PostMerger.py -i Merge_list -o OUTPUT
    ex) python PostMerger.py -i merge_list -o toy_output
    
    ######## Merge_list ########
    <VCF file1 to merge> weight 
    <VCF file2 to merge> weight

    ex)
    /path/to/toy1.vcf 0.5
    /path/to/toy2.vcf 0.5


***


## (3) Citation

TBD.
                              

***


## (4) License
The ***PostMerger*** Software is freely available for non-commercial academic research use. For other usage, one must contact Buhm Han (BH) at buhm.han@snu.ac.kr (patent pending). WE (Hyunjoon Lim, BH) MAKE NO REPRESENTATIONS OR WARRANTIES WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED HERE UNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED "AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE IS UNDERTAKEN AT YOUR OWN RISK. TO THE FULLEST EXTENT ALLOWED BY APPLICABLE LAW, IN NO EVENT SHALL WE BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT, PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT FORESEEABLE AND WHETHER OR NOT WE HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS OR OTHER FINANCIAL LOSS.
