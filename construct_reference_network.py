import sys
import multiprocessing
import numpy as np

param={}
for i in range(1,len(sys.argv)):
    t=sys.argv[i].split("=")
    print(t)
    param[t[0].lower()]=t[1]
    
help_msg="""
usage: python construct_reference_network.py -process=process_value -proportion=proportion_value -PCC=threshold_PCC -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-process: namesumber of work processes used
-proportion: the genes whose expression value was 0 in more than proportion samples/cells
-PCC : set the threshold of PCC [0..1], if the -pvalue set 0, all edges will be outputted to the background_network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the PSSN
-out : the file to store the background_network
"""
if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)

if "-ref" not in param.keys() or "-sample" not in param.keys() or "-process" not in param.keys() or "-out" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()
reference_file=param["-ref"]
sample_file=param["-sample"]
process=int(param["-process"])
file_name=param["-out"]

if "-proportion" in param.keys():
    proportion=float(param["-proportion"])
    if proportion < 0 or proportion > 1:
        print("Please set the correct threshold of p-value in -proportion [0..1]")
        exit()

if "-PCC" in param.keys():
    pcc=float(param["-PCC"])
    if pcc < 0 or pcc > 1:
        print("Please set the correct threshold of p-value in -pvalue [0..1]")
        exit()

def work(i,genes,ref_values):
    print("inner process",i)
    ref_values=np.array(ref_values)
    ref_values=np.corrcoef(ref_values)
    ref_values=np.abs(ref_values)
    for j in range(i+1,len(genes)):
        if ref_values[j][i]<pcc:
            continue
        gene_add=np.add(ref_values[i,0:],ref_values[j,0:])
        gene_add= gene_add.tolist()
        gene_add=list(zip(gene_add,genes))
        gene_add.remove(gene_add[j])
        gene_add.remove(gene_add[i])
        gene_add=sorted(gene_add,key=lambda d:d[0],reverse=True)
        gene=[gene_add[k][1] for k in range(5)]
        gene_partial=" ".join(gene)
        fw=open(file_name,"a")
        fw.write(genes[i] + "_" + genes[j] + "\t" + gene_partial + "\n")
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
    for p in range(len(genes)):
        normal_gene_array=ref[genes[p]]
        tumor_gene_array=data[genes[p]]
        normal_count= normal_gene_array.count(0.0)
        tumor_count=tumor_gene_array.count(0.0)
        if normal_count/len(ref[genes[p]]) > proportion or  tumor_count/len(data[genes[p]]) > proportion:
            gene_pop.append(genes[p])

    for p in range(len(gene_pop)):
        ref.pop(gene_pop[p])
        data.pop(gene_pop[p])
    genes=list(ref.keys())  

    ref_values=list(ref.values())

    pool=multiprocessing.Pool(process)
    
    for i in range(len(genes)-1):

        pool.apply_async(work,args=(i,genes,ref_values,))
    print('Waiting for all subprocesses done...')    
    pool.close()
    pool.join()       

