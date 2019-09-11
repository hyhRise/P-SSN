import os
import sys
import time
import numpy as np
import scipy.stats as stat
import multiprocessing


param={}
for i in range(1,len(sys.argv)):
    t=sys.argv[i].split("=")
    param[t[0].lower()]=t[1]

help_msg="""
usage: python construct_single_network.py -proportion=proportion_value -process=process_value -pvalue=threshold_p_value -background=background_network_file -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-proportion: the genes whose expression value was 0 in more than proportion samples/cells
-process: namesumber of work processes used
-pvalue : set the threshold of p-value [0..1], if the -pvalue set 1, all edges will be outputted to the PSSN
-background : background network to calculate the deltaPCC of edges based on the network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the PSSN
-out : the directory to store the SSN
"""
if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-ref" not in param.keys() or "-sample" not in param.keys() or "-background" not in param.keys() or "-process" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()
reference_file=param["-ref"]
sample_file=param["-sample"]
background=param["-background"]
process=int(param["-process"])

if "-proportion" in param.keys():
    proportion=float(param["-proportion"])
    if proportion < 0 or proportion > 1:
        print("Please set the correct threshold of p-value in -proportion [0..1]")
        exit()

if "-pvalue" in param.keys():
    p_value=float(param["-pvalue"])
    if p_value < 0 or p_value > 1:
        print("Please set the correct threshold of p-value in -pvalue [0..1]")
        exit()

if "-out" not in param.keys():    
    fold="."
else:
    fold=param["-out"]

if not os.path.exists(fold):
    os.mkdir(fold)

def partial_deltapcc(x,y,z,xxx,yyy,zzz):#x,y,z是xxx，yyy，zzz分别是添加的一个元素（就是添加一个样本对应的值）
    xx=np.array(xxx)
    yy=np.array(yyy)
    zz=np.array(zzz)
    xz=stat.linregress(zz,xx)
    yz=stat.linregress(zz,yy)
    rxx=xx-(xz.slope*zz+xz.intercept)
    ryy=yy-(yz.slope*zz+yz.intercept)
    rx=x-(xz.slope*z+xz.intercept)
    ry=y-(yz.slope*z+yz.intercept)
    rxxx=np.append(rxx,rx)
    ryyy=np.append(ryy,ry)
    return (stat.pearsonr(rxxx,ryyy)[0]-stat.pearsonr(rxx,ryy)[0],stat.pearsonr(rxx,ryy)[0])

def ssn_score(deta,pcc,nn):
    if pcc==1:
        pcc=0.99999999
    if pcc==-1:
        pcc=-0.99999999
    z=deta/((1-pcc*pcc)/(nn-1))
    return z 


def parallel_procedure(index,names,gene_filter,gene_correspondence,data,ref,ref_count):
	print("inner process",index)
	fw=open(fold+os.sep+"ppi_"+ names[index]+".txt","w")
	for i in range(len(gene_filter)):
		s=gene_filter[i].split("_")
		pvalue_5=[]
		for j in gene_correspondence[gene_filter[i]]:
			r,r1=partial_deltapcc(data[s[0]][index],data[s[1]][index],data[j][index],ref[s[0]],ref[s[1]],ref[j])
			z=ssn_score(r,r1,ref_count)
			pvalue=1-stat.norm.cdf(abs(z))
			pvalue_5.append(pvalue)
		pvalue_5=np.array(pvalue_5)
		pvalue_threshold=np.array([p_value,p_value,p_value,p_value,p_value])
		compare_result=list(pvalue_5 < pvalue_threshold)
		if compare_result.count(True)==5:
			fw.write(s[0]+"\t"+s[1]+"\n")
	fw.close()

if __name__=="__main__":
    ref={}
    ref_count=0
    f=open(reference_file)
    flag=0
    for p in f:
        flag+=1
        t=p.split()
        if flag==1:
            ref_count=len(t)
            continue
        if t[0] not in ref.keys():
            ref[t[0]]=[float(t[i]) for i in range(1,len(t))]
        else:
            print("Error in ",t[0])
    f.close()
    genes=list(ref.keys())

    data={}
    f=open(sample_file)
    flag=0
    for p in f:
        flag+=1
        t=p.split()
        if flag==1:
            names=t[0:]
            continue
        if t[0] not in data.keys():
            data[t[0]]=[float(t[i]) for i in range(1,len(t))]
        else:
            print("Error in ",t[0])
    f.close()

    gene_pop=[]
    for i in range(len(genes)):
        normal_gene_array=ref[genes[i]]
        tumor_gene_array=data[genes[i]]
        normal_count= normal_gene_array.count(0.0)
        tumor_count=tumor_gene_array.count(0.0)
        if normal_count/len(ref[genes[i]]) > proportion or  tumor_count/len(data[genes[i]]) > proportion:
            gene_pop.append(genes[i])

    for i in range(len(gene_pop)):
        ref.pop(gene_pop[i])
        data.pop(gene_pop[i])
	
    gene_correspondence = {}
    f=open(background)
    flag=0
    for p in f:
        t = p.split()
        if t[0] not in gene_correspondence.keys():
            gene_correspondence[t[0]]=[t[i] for i in range(1,len(t))]
        else:
            print("Error in",t[0])
    f.close()
    gene_filter=list(gene_correspondence.keys())

    pool=multiprocessing.Pool(process)
    for index in range(len(names)):
        pool.apply_async(parallel_procedure,(index,names,gene_filter,gene_correspondence,data,ref,ref_count,))
    print('Waiting for all subprocesses done...')  
    pool.close()
    pool.join()