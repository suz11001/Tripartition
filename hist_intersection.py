import math
import re
import argparse
import os
import numpy as np
import sys
import glob 
import scipy.stats as st
import dendropy
from dendropy import Tree
from dendropy.calculate import treecompare
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')

parser = argparse.ArgumentParser(
     prog='createGFF.py',
     usage='''python createGFF.py --sample [gene sample file bootstrapped directory] --raxml [name of raxml file for sample]''',
     description='''calculate 95% CI for a sample''',
     epilog='''It requires numpy and biopython libraries''')
parser.add_argument('--sample', type=str, help='The name of the breakpoint file', required=False)
parser.add_argument('--raxml', type=str, help='name of the fasta fle', required=False)

args=parser.parse_args()
raxml_tree=args.raxml
raxml_samples_dir=args.sample

def return_intersection(hist_1, hist_2):
    #intsc=np.sum(np.minimum(x,y))
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    return (m-h, m+h)

def calcRF(t1,t2):
    tns = dendropy.TaxonNamespace()
    tree1=Tree.get_from_path(str(t1),"newick",taxon_namespace=tns)
    tree2=Tree.get_from_path(str(t2),"newick",taxon_namespace=tns)
    tree1.encode_bipartitions()
    tree2.encode_bipartitions()
    return (treecompare.unweighted_robinson_foulds_distance(tree1, tree2))

def compareNeighbors(sample):
    samplePath="/home/FCAM/szaman/researchMukul/pgtr/nov8_2020/window_analysis_60perc/02_trees/sequences_control_domain60perc/control/"+str(sample)
    # r=root, d=directories, f = file
    # print(samplePath)                                                                                       

    window1_tree=samplePath+"/window1/RAxML_bestTree.window1"
    window2_tree=samplePath+"/window2/RAxML_bestTree.window2"
    window3_tree=samplePath+"/window3/RAxML_bestTree.window3"
    w1w2=calcRF(window1_tree,window2_tree)
    w1w3=calcRF(window1_tree,window3_tree)
    w2w3=calcRF(window2_tree,window3_tree)

    return(w1w2, w1w3, w2w3)

def generateDistribution(sample,best_tree_rfs):
    print 'sample is: ' + str(sample)
    w1w2s=np.zeros((100))
    w1w3s=np.zeros((100))
    w2w3s=np.zeros((100))

    cwd="/home/FCAM/szaman/researchMukul/pgtr/nov8_2020/window_analysis_60perc/03_bootstrap_trees/sequences_control_domain60perc/control/"

    index=0

    for bootstrap in glob.glob(cwd+str(sample)+'/window1/x*'):
        filename=(os.path.basename(bootstrap))
        w2_tree=cwd+str(sample)+'/window2/'+filename
        w3_tree=cwd+str(sample)+'/window3/'+filename
        w1w2=calcRF(bootstrap,w2_tree)
        w1w3=calcRF(bootstrap,w3_tree)
        w2w3=calcRF(w2_tree,w3_tree)

        w1w2s[index]=w1w2
        w1w3s[index]=w1w3
        w2w3s[index]=w2w3

        index=index+1
       
    maximum=max(max(w1w2s),max(w1w3s),max(w2w3s))+1
    minimum=min(min(w1w2s),min(w1w3s),min(w2w3s))
    nbin=int(maximum-minimum)
    print('number of bins', nbin)

    hist_w1w2s, b1 = np.histogram(w1w2s, bins=nbin, range=(minimum,maximum))
    hist_w1w3s, b2 = np.histogram(w1w3s, bins=nbin, range=(minimum,maximum))
    hist_w2w3s, b3 = np.histogram(w2w3s, bins=nbin, range=(minimum,maximum))

    #create 95% confidence interval for population mean weight                                               
    #print(rfs)
    # print("stats")
    # print(np.mean(w1w2s))
    # print(np.std(w1w2s))
    # print(float(np.std(w1w2s)/math.sqrt(100)))
    # z=(float(np.std(w1w2s)/math.sqrt(100)))
    # print(np.mean(w1w2s)-(z*1.96))
    # print(np.mean(w1w2s)+(z*1.96))
    
    print("w1w2, w1w3, w2w3")
    print(best_tree_rfs)

    print('average of w1w2 bootstraps', w1w2s.mean())
    print('average of w1w3 bootstraps', w1w3s.mean())
    print('average of w2w3 bootstraps', w2w3s.mean())

    print('standard dev of w1w2 bootstraps', np.std(w1w2s))
    print('standard dev of w1w3 bootstraps', np.std(w1w3s))
    print('standard dev of w2w3 bootstraps', np.std(w2w3s))
    
    print("w1w2, w1w3, w2w3")

    x=None

    # if max is w1w3 and min is w1w2
    if  best_tree_rfs.index(max(best_tree_rfs)) == 1 and best_tree_rfs.index(min(best_tree_rfs))==0:
        x=return_intersection(hist_w1w3s,hist_w1w2s)

    # if max is w1w3 and min is w2w3
    elif best_tree_rfs.index(max(best_tree_rfs)) == 1 and best_tree_rfs.index(min(best_tree_rfs))==2:
        x=return_intersection(hist_w1w3s,hist_w2w3s)

    # if max is w2w3 and min is w1w2
    elif  best_tree_rfs.index(max(best_tree_rfs)) == 2 and best_tree_rfs.index(min(best_tree_rfs))==0:
        x=return_intersection(hist_w1w2s,hist_w2w3s)

    # if max is w2w3 and min is w1w3
    elif best_tree_rfs.index(max(best_tree_rfs)) == 2 and best_tree_rfs.index(min(best_tree_rfs))==1:
        x=return_intersection(hist_w1w3s,hist_w2w3s)

    # if max is w1w2 and min is w1w3   
    elif best_tree_rfs.index(max(best_tree_rfs)) == 0 and best_tree_rfs.index(min(best_tree_rfs))==1:
        x=return_intersection(hist_w1w2s,hist_w1w3s)

    # if max is w1w2 and min is w2w3
    elif best_tree_rfs.index(max(best_tree_rfs)) == 0 and best_tree_rfs.index(min(best_tree_rfs))==2:
        x=return_intersection(hist_w2w3s,hist_w1w2s)
    
    if x!=None:
        if x < 0.5:
            print 'sample ' + str(sample) + ' passed the statistical test with ' + str(x)
        else:
            print 'sample ' + str(sample) + ' failed the statistical test with ' + str(x)
    else:
        print('all rfs are equivalent - not possible to distinguish')


    plt.hist(w1w2s, bins=nbin, label='w1w2')
    plt.hist(w1w3s, bins=nbin, label='w1w3')
    plt.hist(w2w3s, bins=nbin, label='w2w3')
    plt.legend()
    plt.savefig('hist_Transfer_'+str(sample)+'.png')

def loopsample(transfers_dict):
    transfers=transfers_dict
    for i in range(1,101):                                                                                                                                                                              
        print(i)                                                                                                                                                                                        
        if transfers.get(i) != 0:                                                                                                                                                                       
            print('this sample '+ str(i) +' contains domain transfer events')                                                                                                                          
            x=compareNeighbors(i)                                                                                                                                                                       
            generateDistribution(i,x)                                                                                                                                                                   
        else:                                                                                                                                                                                           
            print('this sample '+ str(i) +' does not contain any domain transfer events')                                                                                                               
            x=compareNeighbors(i)                                                                                                                                                                       
            generateDistribution(i,x) 

if __name__ == "__main__":
    
    #x=compareNeighbors(11)
    #generateDistribution(11,x)
    domain_samples="/home/FCAM/szaman/researchMukul/pgtr/nov8_2020/domain_transfers.txt"
    transfers = { int(line.split(":")[0]) : int(line.split(":")[1]) for line in open(domain_samples) }
    loopsample(transfers)
    
