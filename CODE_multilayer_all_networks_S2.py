#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 22:10:51 2021

@author: giriraj
"""

'''
we gather all the data and create the layers for multulayer network for further analysis
'''

import os
import pandas as pd
import networkx as nx
import pickle as pk
import backboning as bck
from cdlib import algorithms, ensemble, evaluation, NodeClustering
import ClusterEnsembles as CE
import time
import numpy as np

os.chdir('~/paper2_codes') # set current directory here

cwd = os.getcwd()

'''
First you prune the coexpression network and then create the subnetworks for PPI, TF and ncRNA networks for the same
This is being done only for SOD1 and TDp43 as they are only ones which will be analyzed later.

'''

def backboner(gene_file):
    '''gene_file: either file name or dataframe'''
    
    if type(gene_file) == str:
        gene_net = pd.read_csv(cwd+'/datasets/' + gene_file , sep=' ',names=['source','target','weight'])
        b_net = gene_net[['source','target','weight']]
        b_net = b_net[pd.to_numeric(b_net['weight'], errors='coerce').notnull()] #to check any non-numeric type weights
        
        try: b_net['weight'] = b_net['weight'].apply(lambda x: float(x))
        except: pass
        
        b_net = b_net[b_net['weight'] > 0.6] #filtering edges with RLowPC weight < 0.6
        b_net.columns = ['src','trg','nij']
        b_net_table = bck.noise_corrected(b_net)
    else:
        b_net = gene_file[['source','target','weight']]
        b_net = b_net[pd.to_numeric(b_net['weight'], errors='coerce').notnull()] #to check any non-numeric type weights
        
        try: b_net['weight'] = b_net['weight'].apply(lambda x: float(x))
        except: pass
        
        b_net = b_net[b_net['weight'] > 0.6] #filtering edges with RLowPC weight < 0.6
        b_net.columns = ['src','trg','nij']
        b_net_table = bck.noise_corrected(b_net)
    
    val = 0; thr = 1
    while val == 0:
        threshold = b_net_table['score'].mean() + thr*b_net_table['score'].std()
        b_net_backbone = bck.thresholding(b_net_table, threshold)
        chk_graph = nx.from_pandas_edgelist(b_net_backbone,source='src',target='trg',edge_attr='nij')
        
        if nx.is_connected(chk_graph):
            thr += 0.5
        elif not nx.is_connected(chk_graph):
            thr -= 0.5
            threshold = b_net_table['score'].mean() + thr*b_net_table['score'].std()
            b_net_backbone = bck.thresholding(b_net_table, threshold)
            val = 1

    gene_df = set(b_net_backbone['src'].values.tolist() + b_net_backbone['trg'].values.tolist())    
    return b_net_backbone, gene_df



def net_create(df,net):
    try:
        x = nx.subgraph(net,df); y = max(nx.connected_components(x), key=len)
        x = nx.subgraph(x,y)
        x = nx.to_pandas_edgelist(x,source='source',target='target')
        return x
    except: print('Network not created :( ')



st1 = time.time()


def node_dic_maker(df,src='source',trg='target'):
    a = [tuple(i) for i in df[[src,'srcInd']].values.tolist()]
    b = [tuple(i) for i in df[[trg,'trgInd']].values.tolist()]
    return dict(a + b)
    

def table_maker(dat):
    '''generates a beautiful table of the community optimization performed before'''
    from prettytable import PrettyTable
    import numpy as np
    tabx = PrettyTable(); tmp = []
    for item in dat:
        for k,v in item[1].items():
            var = v[1].communities
            p = len(var)
            var = [len(t) for t in var]; q = round(np.mean(var),4); r = round(np.std(var),4)
            s = round(v[0].score,4)
            tmp.append([item[0],k,p,q,r,s])
    tabx.field_names = ['Layer','Algo','Total COMS','COM.Size(Avg.)', 'COM.Size(Std.)','Mod.Density.Score']
    tabx.add_rows(tmp)
    
    return tabx


def consensus_maker(comList,nodeDict=None,lname=None,comGraph=None,multi=False,partition_size=None):
    '''code to identify the consensus clustering among the different layers of the network.
    This will take  the inputgenerated from the CDLIB library and optimized for the maximization of 
    modularity density
    -- comList :: list of NodeClustering Objects, same net partitioned multiple times or multilayer partitions
    -- nodeDict :: dictionary of node names and unique indices
    -- lname :: name of the network being processed for identification
    -- comGraph :: graph associated with the network, for multilayer its nx.compose() output
    -- multi: Boolean, checks whether you are running it for multilayer creation or not
    '''
    
    if not multi:
        import operator
    
        mat = np.empty((len(comList),len(nodeDict))); mat[:] = np.nan
        comList = [x.communities for x in comList]
        
        for i in range(len(comList)):
            values = comList[i]; keys = list(range(len(values)))
            for k,v in zip(keys,values):
                for idx in v: mat[i][idx] = k #for idx to work the node indices must be numeric
        
        #consensus partions are generated for range between the size of smallest community to the  largest community
        consensus_clusters = {}; consensus_scores = {}
        uv = len(max(comList,key=len)) #getting min and max sized communities in list of communities
    
        if uv !=1:
            tmp = CE.cluster_ensembles(mat,verbose = False, nclass = uv, solver = 'all',random_state=108)
            #formatting the cluster array to convert it into list of lists
            vx = [[] for p in range(max(tmp)+1)] #empty list of list , length same as tmp; CE returns a list of community assignments to  all nodes
            for p in range(len(tmp)): vx[tmp[p]].append(p)
            consensus_clusters[int(uv)] = NodeClustering(vx, graph = comGraph, method_name = lname)
            consensus_scores[int(uv)] = evaluation.modularity_density(comGraph, consensus_clusters[int(uv)])
        
        else: #when there is only one cluster and no partitions found by the CDLIB algorithm
            tmp = CE.cluster_ensembles(mat,verbose = False, nclass = 2, solver = 'all',random_state=108)
            vx = [[] for p in range(max(tmp)+1)] #empty list of list , length same as tmp
            for p in range(len(tmp)): vx[tmp[p]].append(p)
            consensus_clusters[0] = NodeClustering(vx, graph = comGraph, method_name = lname)
            consensus_scores[0] = evaluation.modularity_density(comGraph, consensus_clusters[0])
            
            
        idx = max(consensus_scores.items(), key=operator.itemgetter(1))[0]
        return consensus_clusters[idx] , consensus_scores[idx]



    if multi:
        '''in this case it will take comList with node names and format is different'''
        try:
            print('Running consensus in multilayer mode')
            
            x1 = []
            for i in comList:
                x1.append(sum(i,[])) # summing list of lists with and empty list just joins up the list, abstract datatype
            
            x1 = list(set(sum(x1,[]))); vals = list(range(len(x1)))
            md = dict(zip(x1,vals))
            
            mat = np.empty((len(comList),len(md))); mat[:] = np.nan
            
            for i in range(len(comList)): #layer level
                for j in range(len(comList[i])):
                    for k in comList[i][j]:
                        idx = md[k]
                        mat[i][idx] = j
    
            consensus_clusters = {}
            
            if partition_size:
                for j in range(min(partition_size) , max(partition_size)): #run consensus for min to max partitions sizes
                    tmp = CE.cluster_ensembles(mat,verbose = False, nclass = j, solver = 'all',random_state=108)
                    #formatting the cluster array to convert it into list of lists
                    vx = [[] for p in range(max(tmp)+1)] #empty list of list , length same as tmp; CE returns a list of community assignments to  all nodes
                    for p in range(len(tmp)): vx[tmp[p]].append(p)
                    consensus_clusters[j] = vx
                    
            if not partition_size:
                print('enter partition sizes of all layers for multilayer mode; it should be a list of sizes')
                raise SystemExit
                    
            return consensus_clusters, md
        except Exception as e:
            print(e)
            print('Please check the format of the community List passed, it should not be numeric')
        
        
        


def oslom_process(layerDF,src='source',trg='target',wgt='weight'):    
    from argparse import Namespace
    import oslom
    import networkx as nx
    
    args = Namespace()
    args.min_cluster_size = 0
    args.oslom_exec = oslom.DEF_OSLOM_EXEC
    args.oslom_args = oslom.DEF_OSLOM_ARGS
    
    #snippet to check the weight and column names as needed
    try: edge = [tuple(i) for i in layerDF[[src,trg,wgt]].values.tolist()]
    except:
        layerDF['weight'] = 1 # setting weight for oslom run 1 for all edges, will see later if actual weight to be added
        edge = [tuple(i) for i in layerDF[[src,trg,'weight']].values.tolist()]
    
    clusters = oslom.run_in_memory(args,edge)
    
    #organising results
    nd = node_dic_maker(layerDF,src=src,trg=trg)
    var = clusters[0]['clusters']
    clus = []
    for p in var:
        tmp = []
        for q in p['nodes']:
            tmp.append(nd[q['id']])
        clus.append(tmp); tmp = []
    
    g = nx.from_pandas_edgelist(layerDF,source = src,target = trg,create_using=nx.Graph())
    coms = NodeClustering(clus, graph = g, method_name="oslom")
    
    return coms

def partition_creator(vector):
    
    graphX = vector[0]; net_name = vector[1]; layerDF = vector[2]    
    nd = node_dic_maker(layerDF,layerDF.columns[0],layerDF.columns[1])
    
    st = time.time()
    # INFOMap ALGORithm -----------------------------------------------------------------------------------
    print('InfoMAP Clustering for ',net_name, ' started')
    #we run the algo for 20 times with same parameters as these algos are stochastic and not deterministic
    imap_var = [algorithms.infomap(graphX) for rx in range(20)]
    
    
    com_imap, sco_imap = consensus_maker(imap_var, nd, net_name, graphX)
    
    print('InfoMAP Clustering for ',net_name, ' finished in', round(time.time()-st,2))
    
    
    #OSLOM Algorithm ---------------------------------------------------------------------------------------
    st = time.time()
    print('OSLOM Clustering for ',net_name, ' started')
    #we run the algo for 20 times with same parameters as these algos are stochastic and not deterministic
    cols = layerDF.columns
    
    if 'source' and 'target' in cols:
        oslom_var = [oslom_process(layerDF) for rx in range(20)]
        com_oslom, sco_oslom = consensus_maker(oslom_var, nd, net_name, graphX)
        
    elif 'src' and 'trg' and 'score' in cols:
        oslom_var = [oslom_process(layerDF,src='src',trg='trg',wgt='score') for rx in range(20)]
        com_oslom, sco_oslom = consensus_maker(oslom_var, nd, net_name, graphX)
    print('OSLOM Clustering for ',net_name, ' finished in', round(time.time()-st,2))
    
    
    #DER clustering (supports weighted network)--------------------------------------------------------------------------------------
    print('DER Clustering for ',net_name, ' started')
    st = time.time()
    walk_len = ensemble.Parameter(name="walk_len",start=3,end=12,step=2)
    com_DER, sco_DER = ensemble.grid_search(graph=graphX,
                                            method=algorithms.der,
                                            parameters=[walk_len],
                                            quality_score=evaluation.modularity_density,
                                            aggregate=max)
    
    print('DER Clustering for ',net_name, ' finished in ', round(time.time()-st,2))
    
    
    #MARKOV clustering -----------------------------------------------------------------------------------
    print('Markov Clustering for ',net_name, ' started')
    st = time.time()
    inflation = ensemble.Parameter(name="inflation", start=1.2, end=2.2, step=0.4)
    expansion = ensemble.Parameter(name="expansion", start=2, end=5, step=1)
    com_markov, sco_markov = ensemble.grid_search(graph=graphX, method=algorithms.markov_clustering,
                                                  parameters=[inflation,expansion],
                                                  quality_score=evaluation.modularity_density,
                                                  aggregate=max)
    
    
    print('Markov Clustering for ',net_name, ' finished in ', round(time.time()-st,2))

    #louvain clustering -----------------------------------------------------------------------------------
    print('Louvain Clustering for ',net_name, ' started')
    resolution = ensemble.Parameter(name="resolution", start=0.1, end=1, step=0.2)
    randomize = ensemble.BoolParameter(name="randomize")
    com_louv, sco_louv = ensemble.grid_search(graph=graphX, method=algorithms.louvain,
                                                     parameters=[resolution, randomize],
                                                     quality_score=evaluation.modularity_density,
                                                     aggregate=max)
    
    print('LOUVAIN Clustering for ',net_name, ' finished in', round(time.time()-st,2))


    #RB_POTS ----------------------------------------------------------------------------------------------
    st = time.time()
    print('RB_POTS Clustering for ',net_name, ' started')
    res = ensemble.Parameter(name='resolution_parameter',start=1, end = 10, step = 1)
    com_pot, sco_pot = ensemble.grid_search(graph=graphX,method=algorithms.rb_pots,
                                            parameters=[res],
                                            quality_score=evaluation.modularity_density,
                                            aggregate=max)
    
    print('RB_POTS Clustering for ',net_name, ' finished in ', round(time.time()-st,2))
    
    
    # ------------------------------------------------------------------------------------------------------
    
    com_set = {'MARKOV':[sco_markov, com_markov],'LOUVAIN':[sco_louv,com_louv],
               'IMAP':[sco_imap, com_imap],'DER':[sco_DER, com_DER],
               'RB_POTS':[sco_pot, com_pot],'OSLOM':[sco_oslom, com_oslom]}
    

    return net_name, com_set, nd

#funtion to assign numeric type nodes as some of the algos require this, so making it universal for all 
def node2numeric(df,source,target):
    nodes = list(set(df[source].values.tolist()+ df[target].values.tolist()))
    values = list(range(len(nodes)))
    df_dict = dict(zip(nodes,values))
    df['srcInd'] = df[source].map(df_dict); df['trgInd'] = df[target].map(df_dict)
    return df





#%%

'''start from here loading the main networks'''

'''--------------------------- creating individual multilayer network for every dataset --------------------------------------'''


'''
This part of the code will run the netowrk clustering using CDLIB as ensemble run to get the best partitions.
Grid Search methodology will be used
'''

#load previous network data first
   
with open(cwd+'/datasets/begin_nets','rb') as fp:
    rna_net,tf_net,ppi_net = pk.load(fp)
 
#NOTE: Comment every name except the file id to be processed.
#gene_file = 'GSE158264_coex_net.csv' #should be string format, will be used to save final community file name
#gene_file = 'tdp43_coex_net.csv'
#gene_file = 'sod_coex_net.csv'
#gene_file = 'GSE145677_coex_net.csv'
gene_file = 'GSE40438_coex_net.csv'




ppi_network = nx.from_pandas_edgelist(ppi_net,source='source',target='target')

tf_network = nx.from_pandas_edgelist(tf_net,source='source',target='target',edge_attr='weight')

rna_network = nx.from_pandas_edgelist(rna_net,source='source',target='target',edge_attr='weight')



gene_coex_net, gene_df = backboner(gene_file) #pruning network

gene_ppi_net = net_create(gene_df,ppi_network)

gene_tf_net = net_create(gene_df,tf_network)
gene_tf_net, tf_df = backboner(gene_tf_net) #pruning network

gene_rna_net = net_create(gene_df,rna_network)
gene_rna_net, rna_df = backboner(gene_rna_net) #pruning network



#%%

'''
This part of the code reads the subnetworks generated using random walk algo in R and then generate a ALS specific gene 
neighborhood fold change plot using the PRISM package
'''

# import glob


# files = glob.glob(cwd+'datasets/DEGs/networkData26Jan2021/coex_subnet_SOD1/*')
# main_df = pd.read_csv(cwd+'datasets/DEGs/modified_sod_expData.csv')
# geneDF = {}; geneDF2 = {}; finalDF = pd.DataFrame()

# for f in files:
#     key = f.split('_')[-1] #taking out the gene names from the file names
    
#     df = pd.read_csv(f)
#     df = df.sort_values(['Score'],ascending=False).head(30) #selecting top 30 by score, as done in the original paper of RW method
    
#     df_main = main_df[['Symbol','log2FoldChange','pvalue','padj']]
#     value = pd.merge(df,df_main,left_on='NodeNames',right_on='Symbol').drop(labels='NodeNames',axis=1).reset_index(drop=True)
    
#     geneDF[key] = value
    
#     finalDF[key] = value['log2FoldChange'].values
    
#     finalDF.to_csv(cwd+'datasets/DEGs/networkData26Jan2021/subnet_mapped_sod1/subnet_foldChange_sod1.csv',index=False)
        

#%%


#assigning numeric indices to the nodes of the graphs
gene_coex_net = node2numeric(gene_coex_net,'src','trg')
gene_ppi_net = node2numeric(gene_ppi_net,'source','target')
gene_rna_net = node2numeric(gene_rna_net,'src','trg')
gene_tf_net = node2numeric(gene_tf_net,'src','trg')


#generating graphs
gene_coex_net1 = nx.from_pandas_edgelist(gene_coex_net,source='srcInd',target='trgInd',create_using=nx.Graph,edge_attr='nij')
gene_ppi_net1 = nx.from_pandas_edgelist(gene_ppi_net,source='srcInd',target='trgInd',create_using=nx.Graph)
gene_tf_net1 = nx.from_pandas_edgelist(gene_tf_net,source='srcInd',target='trgInd',create_using=nx.Graph,edge_attr='nij')
gene_rna_net1 = nx.from_pandas_edgelist(gene_rna_net,source='srcInd',target='trgInd',create_using=nx.Graph,edge_attr='nij')


in_nets2 = [ (gene_ppi_net1,'ppi',gene_ppi_net), (gene_coex_net1,'coex',gene_coex_net),
           (gene_tf_net1,'tf',gene_tf_net) , (gene_rna_net1,'rna',gene_rna_net)]


dat_gene1 = []
for inpt in in_nets2:
    dat_gene1.append(partition_creator(inpt))

#%%
'''
code to identify the best of the 5 community detection algorithms for each layer of network; and rank the communities
'''
best_partitions = {}; partition_size = []

for i in dat_gene1:
    tmp = 0; val = ()
    for k,v in i[1].items():
        if v[0].score > tmp:
            tmp = v[0].score; val = (k , v[0].score , i[2] , v[1].communities)
    best_partitions[i[0]] = val
    partition_size.append(len(val[3]))


'''
The best cluster of the consensus partitions of different cluster sizes is evaluated using the modularity score
The cluster/class size which gives the maximum modularity is choosen, if tie, then the partition with higher partitions is choosen
'''
'''
based on the algorithm comparisons we identify the best algo for each case to do the community partition
generating the consensus in the all 4 layer partitions data. The consensus is optimized against the
aggregate network generated using the nx.compose method
consensus_maker(comList,nodeDict,lname,comGraph)
'''

gene1_finalP = {}

for ky,va in best_partitions.items():
    gene1_finalP[ky] = [va[2],va[3]]


#remaps the node indices to the node names

for k,v in gene1_finalP.items():
    vals = list(v[0].keys()); keys = list(v[0].values())
    
    ndc = dict(zip(keys,vals))
    
    for j in range(len(v[1])):
        tmp = [ndc[x] for x in v[1][j]]
        v[1][j] = tmp
    
    gene1_finalP[k] = v[1]
    

'''
Evaluating the consensus matrix of the node pairs
Steps:
    First make the community matrix for every layer, based on whether a given node pair exist in the same community for a layer
    collect all nodes and create a dictionary
    then loop over the unique nodes for each layer partition
    if the node pair exists within any partition then its 1, else 0
    update the individual community matrix
    take average of all 4 community matrices which gives the consensus matrix of the gene type, tdp43 and sod1
'''
def community_matrix(nested_gene_list,layerKey):
    '''
    nested_gene_list: dict of layers as keys and partitioned nested list of node names
    layerKey: one of the 'coex', 'tf', 'rna' and 'ppi'
    '''
    from itertools import chain
    import numpy as np
    
    #collect all genes to create dict
    glist = chain(*list(chain(*list(nested_gene_list.values())))) #note: two chain function required. inner nested loop and dict nested loop
    glist = sorted(list(set(list(glist))))
    
    matx = np.zeros((len(glist),len(glist)))
    dicx = dict(zip(glist,list(range(len(glist)))))
    
    for i in glist:
        for j in nested_gene_list[layerKey]:
            if i in j:
                for k in j:
                    if k != i: #taking care of diagonals which is obsolete
                        matx[dicx[i],dicx[k]] = 1    
    return matx, dicx


klist = ['coex','tf','rna','ppi']
smat = []; tmat = [];

for key in klist:
    tmp2,td2 = community_matrix(gene1_finalP,key)
    smat.append(tmp2)

#average of all four layers is the consensus matrix probability of each pair of nodes
smat = sum(smat)/4.0
tmat = sum(tmat)/4.0


# '''

# here the multilayers are collection of 4 different layers and thus nx.compose graph doesn't reflect a proper
# graph to evaluate the modularity density metric. some other metric is to be used which doesn't require graph'

# TO BE MODIFIED

# '''

y = [v for k,v in gene1_finalP.items()]

multi_net_gene = consensus_maker(comList = y ,multi=True, partition_size= partition_size)

#arranging the partition assignments

def updater2(datt):
    ky = list(datt[1].values()); vl = list(datt[1].keys())
    datt[1] = dict(zip(ky,vl))    
    for k,v in datt[0].items():
        for i in range(len(v)):
            v[i] = [datt[1][j] for j in v[i]]
        datt[0][k] = v
    return datt


multi_net_gene = updater2(list(multi_net_gene))

#generating community matrix for all cluster sizes of multilayer consensus partition
con_gene = {}
for key,val in multi_net_gene[0].items():
    tmp,mtd1 = community_matrix(multi_net_gene[0],key)
    con_gene[key] = tmp


#evaluating similarity between the community matrix of multi-consensus partitions and the consensus matrix of individual layers
# the scores will be considered significant upto 2 decimal places only. By this logic the min frobenius distance score


similarity_mat = {}
for key,val in con_gene.items():
    tmp = round(abs(np.linalg.norm(smat,ord='fro')-np.linalg.norm(val,ord='fro')),4)
    similarity_mat[key] = tmp
    

'''
saving the best consensus and coexpression partitions
'''

consIDX = min(similarity_mat , key = similarity_mat.get)
coex_best = best_partitions['coex'][3]

#writing best coex default partitions
with open(cwd+'/datasets/'+gene_file.split('.')[0]+'_coex_partitions.csv','w') as fp:
    for i in coex_best:
        for j in i:
            fp.write(j+',')
        fp.write('\n')

#writing best coex consensus partitions
with open(cwd+'/datasets/'+gene_file.split('.')[0]+'_selected_partition.csv','w') as fp:
    for i in multi_net_gene[0][consIDX]:
        for j in i:
            fp.write(j+',')
        fp.write('\n')
