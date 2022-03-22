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

def return_intersection(hist_1, hist_2):
    #intsc=np.sum(np.minimum(x,y))
    minima = np.minimum(hist_1, hist_2)
    intersection = np.true_divide(np.sum(minima), np.sum(hist_2))
    return intersection

def calcRF(t1,t2):
    tns = dendropy.TaxonNamespace()
    tree1=Tree.get_from_path(str(t1),"newick",taxon_namespace=tns)
    tree2=Tree.get_from_path(str(t2),"newick",taxon_namespace=tns)
    tree1.encode_bipartitions()
    tree2.encode_bipartitions()
    return (treecompare.unweighted_robinson_foulds_distance(tree1, tree2))

def compareNeighbors(ml_trees_path):
    samplePath=ml_trees_path
    # r=root, d=directories, f = file
    # print(samplePath)                                                                                       

    window1_tree=samplePath+"/RAxML_bestTree.window1"
    window2_tree=samplePath+"/RAxML_bestTree.window2"
    window3_tree=samplePath+"/RAxML_bestTree.window3"
    w1w2=calcRF(window1_tree,window2_tree)
    w1w3=calcRF(window1_tree,window3_tree)
    w2w3=calcRF(window2_tree,window3_tree)

    return(w1w2, w1w3, w2w3)

def generateDistribution(threshold,sample,best_tree_rfs):
    print 'sample is: ' + str(sample)
    w1w2s=np.zeros((100))
    w1w3s=np.zeros((100))
    w2w3s=np.zeros((100))

    index=0

    for bootstrap in glob.glob(str(sample)+'/x*.window1_bs'):
        filename=(os.path.basename(bootstrap)).split(".")[0]
        w2_tree=str(sample)+'/'+filename+".window2_bs"
        w3_tree=str(sample)+'/'+filename+".window3_bs"
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
        if x < threshold:
            print 'sample ' + str(sample) + ' passed the statistical test with ' + str(x)
        else:
            print 'sample ' + str(sample) + ' failed the statistical test with ' + str(x)
    else:
        print('all rfs are equivalent - not possible to distinguish')

def plot(w1w2s,w1w3s,w2w3s):

    plt.hist(w1w2s, bins=nbin, label='w1w2')
    plt.hist(w1w3s, bins=nbin, label='w1w3')
    plt.hist(w2w3s, bins=nbin, label='w2w3')
    plt.legend()
    plt.savefig('histogram.png')

def loopsample(threshold,bootstrap_path, ml_trees_path):
     x=compareNeighbors(ml_trees_path)
     generateDistribution(threshold,bootstrap_path,x)                                                                                                                                                                   


if __name__ == "__main__":
    
     parser = argparse.ArgumentParser(
     prog='hist_intersection.py',
     usage='''python hist_intersection.py --thresh [required overlap for gene families ] --bs_sample [gene sample file bootstrapped directory] --ml_sample [name of raxml file for sample]''',
     description='''determine whether gene family has sub-gene transfer (presence/absence)''',
     epilog='''It requires numpy, dendropy, scipy libraries''')
     parser.add_argument('--thresh', type=str, help='minimum overlap required between distributions', required=True)
     parser.add_argument('--bs_sample', type=str, help='path to bootstrap samples', required=True)
     parser.add_argument('--ml_sample', type=str, help='path to maximum likelihood trees', required=True)
     parser.add_argument('--plot', type=bool, help='plot the RF-score distribition between partitions', required=False)

     args=parser.parse_args()
     threshold=float(args.thresh)
     ml_trees_path=args.ml_sample
     bootstrap_path=args.bs_sample
    
     loopsample(threshold,bootstrap_path, ml_trees_path)
    
