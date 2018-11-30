#注意更改路径
import gc, argparse, sys, os, errno
%pylab inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_style('whitegrid')
import h5py
import os
from tqdm import tqdm_notebook as tqdm
import scipy
import sklearn
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')
 
transcript_table = pd.read_table('/BioII/chenxupeng/student/data/other_annotations/transcript_anno.txt')
 
lengths = (transcript_table.loc[:, 'end'] - transcript_table.loc[:, 'start']).astype('str')
feature_labels = transcript_table.loc[:, 'transcript_id'] + '|' +  transcript_table.loc[:, 'transcript_name'] + '|' + lengths
feature_labels.index = transcript_table.loc[:, 'transcript_id']
feature_labels = feature_labels[~feature_labels.index.duplicated()]
 
mx = pd.read_table('<your path to>/proj_exRNA/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx')
#,index_col=0
mx.head()
 
mx['geneID'] = feature_labels.loc[mx['geneID'].values].values
mx=mx.set_index('geneID')
 
#你需要改变路径到自己文件夹下,scirep_sequential_qc.txt文件在路径/BioII/chenxupeng/student/stat下
ref2 = pd.read_table('/home/chenxupeng/projects/training/data/expression_matrix/GSE71008.txt',sep='\t',index_col=0)
 
samplename = np.array(['Sample_N1','Sample_N7','Sample_N13','Sample_N19','Sample_N25'])
 
ref2.loc[:,samplename].head()
 
mx_nonzero = mx.loc[:,samplename].iloc[np.where(mx.loc[:,samplename].sum(axis=1)!=0)]
#mx_nonzero.head(50000)
 
mx.loc[:,samplename]
 
ref2_nonzero = ref2.loc[:,samplename].iloc[np.where(ref2.loc[:,samplename].sum(axis=1)!=0)]
ref2_nonzero.head(50000)
 
index_overlap = np.intersect1d(ref2_nonzero.index,mx_nonzero.index)
 
#算一下两种方法的index的overlap的占比
index_overlap.shape[0]/mx_nonzero.shape[0], index_overlap.shape[0]/ref2_nonzero.shape[0]  
 
#计算相关系数
pearsonr(np.array(mx_nonzero.loc[index_overlap]).ravel(), np.array(ref2_nonzero.loc[index_overlap]).ravel())
