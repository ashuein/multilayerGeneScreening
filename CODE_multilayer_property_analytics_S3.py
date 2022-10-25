#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 20:22:28 2022

@author: giriraj


Community ranking based on original coexpression network

Enrichr API link

GO terms mapping to the communities

Inter and Intra community impact of control nodes

"""

import json, requests, os
import pandas as pd
import matplotlib.pyplot as plt
from math import gcd



os.chdir('~/paper2_codes') #set current directory here

cwd = os.getcwd()

def gene_enrichment_analysis(glist):
    '''fun to get the enriched pathways in the given gene list'''
    
    assert type(glist) == list, 'you must provide a list of genes (list datatype)'
    #generate Enrichr job process ID -------------------------------------------------
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    
    #this has to be replaced with the genes in your list------------------------------------------>>>>>>>
    genes_str = '\n'.join(glist)
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    processID = json.loads(response.text)['userListId']
    
    #retrieve gene set enrichment analysis terms ------------------------------------
    databases = ['BioPlanet_2019', 'KEGG_2021_Human','Panther_2016','Reactome_2016']
    
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    data = []
    
    for library in databases:
        response = requests.get(ENRICHR_URL + query_string % (processID,library))
        if not response.ok:
            raise Exception('Error fetching enrichment results from '+library)
        data.append(json.loads(response.text))
        print('Library Processed:',library)
        
    #process the retrieved data-------------------------------------------------------
    tmp = []
    for base in data:
        for k,v in base.items():
            for paths in v:
                tmp.append([k,paths[1],paths[2],paths[6],len(paths[5]),paths[5]])
    tmp = list(zip(*tmp))
    try:
        pathwayDF = pd.DataFrame({'Database':tmp[0],'Pathway':tmp[1],'p-val':tmp[2],
                                  'p-adj':tmp[3],'genes_mapped':tmp[4],'genes':tmp[5]})
        
        pathwayDF = pathwayDF.sort_values(by=['p-adj'], ascending=True).reset_index(drop=True)
        return pathwayDF,data
    except:
        print('No enrichment found')
        return None,None
    

def factorization(n):
    factors = []
    def get_factor(n):
        x_fixed = 2
        cycle_size = 2
        x = 2
        factor = 1
        while factor == 1:
            for count in range(cycle_size):
                if factor > 1: break
                x = (x * x + 1) % n
                factor = gcd(x - x_fixed, n)
            cycle_size *= 2
            x_fixed = x
        return factor
    while n > 1:
        next = get_factor(n)
        factors.append(next)
        n //= next
    return factors

#%%

'''read main file here, comment every file name except the one to be processed'''
#fname = 'GSE158264_coex_net'
fname = 'tdp43_coex_net'
#fname = 'sod_coex_net'
#fname = 'GSE145677_coex_net'
#fname = 'GSE40438_coex_net'

#generating dataframe for selected partition analysis
com_file = cwd+'/datasets/'+fname+'_selected_partition.csv'
df = pd.read_csv(com_file,sep='\t',header=None)

#coexpression network dataframe default partitions
com_file = cwd+'/datasets/'+fname+'_coex_partitions.csv'
df2 = pd.read_csv(com_file,sep='\t',header=None)




#%%

'''
This part of the code will parse the set of genes of individual communities to Enrichr API
and then collect the best mapped pathway listing for the genes among the several types of databases KEGG, BioPlanet2019, etc.
convert this into the type required to do the downstream analysis.
'''

#main file consensus
com_paths = {}
for i in range(len(df.index)):
    com_paths['com '+str(i+1)] = gene_enrichment_analysis(df.iloc[i,0].strip(',').split(','))
delk = []
for k,v in com_paths.items(): #removing 'none' items from the dict, due to no enrichment mapping. This will reduce no. of coms counted for plotting later
    try:
        if not v[0]: delk.append(k)
    except: pass
for k in delk: com_paths.pop(k, None)

#main file coexpression
com_paths_coex = {}
for i in range(len(df2.index)):
    com_paths_coex['com '+str(i+1)] = gene_enrichment_analysis(df2.iloc[i,0].strip(',').split(','))
delk = []
for k,v in com_paths_coex.items():
    try:
        if not v[0]: delk.append(k)
    except: pass
for k in delk: com_paths_coex.pop(k, None)

#%%

'''PLOT 1.'''

'''
this part of the code takes the identified enriched pathways in the communities of the consensus network, 
gets the top10 pathways identified, sorts them based on number of genes mapped, finally makes the plot
of pathways mapped, and p-value while color of the dot denoting number of genes mapped
'''

#Note: this code identifies pathways to which ALS nodes have been mapped to

alsOD = pd.read_excel(cwd+'/datasets/als_od.xlsx')

#first read all community groups and filter the GO pathways
go_paths = []; als_map = {}

for k,v in com_paths.items():
    vals = v[0].sort_values(by='p-adj')
    vals = vals.head(10).reset_index(drop=True)
    for i in range(len(vals.loc[:,'genes'])):
        tmp = set(vals.loc[:,'genes'][i]).intersection(set(alsOD['node'].values.tolist())) #identified ALS nodes
        if tmp and vals.loc[:,'p-adj'][i] <= 0.1: #filtering only significant among the mapped als nodes with pathway q-val <=0.1
            go_paths.append([k,vals.loc[:,'Pathway'][i],
                        vals.loc[:,'p-adj'][i],vals.loc[:,'genes_mapped'][i],tmp]) #storing als node community number, identified als nodes,pathway, q-value
   
go_paths = list(zip(*go_paths))

plt.rcParams.update({'font.size': 16})
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["font.family"] = "serif"

coms = go_paths[0]

y = go_paths[1]
x = go_paths[2]
z = [int(i) for i in go_paths[3]]

pcm = plt.scatter(x, y, lw=0.5,c=z, cmap=plt.cm.get_cmap('Dark2', 6),s=300)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.xlabel('adjusted p-value',fontdict=dict(weight='bold'))
plt.title(fname.split('_')[0])

plt.rcParams.update({'font.size': 12})
cbar = plt.colorbar(pcm,ticks=list(range(int(min(z)),int(max(z)+1),int((max(z)+1)/4))),label = 'gene overlaps')
cbar.set_label(label='gene counts',weight='bold')

for i in range(len(x)): plt.annotate(coms[i],(x[i],y[i]),ha='center',va='bottom')


#%%

'''PLOT 2.'''

# This part of code identifies the pathways which are in top 5 in the communities which have been identified to be associated with als nodes
coms = sorted(list(set(coms)))
selected_coms = {}
for i in coms:
    selected_coms[i] = com_paths[i][0].sort_values(by='p-adj').head(5).reset_index(drop=True)
#--------------------------------------------------------------------
#deciding the number of rows and columns for subplots for als associated community top 5 pathways. factorize total plots number and consider three cases of factors
#Case1: no factors: rows = no. of community, cols = 1
#Case2: 2 factors: rows = max(2factors); cols = min(2factors)
#Case3: >2factors: cols = smallest no., rows = product of all numbers except smallest ( 1 is not considered)

chk = factorization(len(coms))
if len(chk) == 1: rows = chk[0]; cols = 1
if len(chk) == 2: rows = max(chk); cols = min(chk)
if len(chk) > 2:
    from numpy import prod
    cols = min(chk); chk.pop(0)
    rows = prod(chk)
#--------------------------------------------------------------------
plt.rcParams.update({'font.size': 14})
fig, axs = plt.subplots(rows,cols, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)
axs = axs.ravel()
for i in range(len(coms)):
    
    y = selected_coms[coms[i]]['Pathway'].values.tolist()
    x = selected_coms[coms[i]]['p-adj'].values.tolist()
    z = [int(i) for i in selected_coms[coms[i]]['genes_mapped'].values.tolist()]
    
    pcm = axs[i].scatter(x, y, lw=0.5,c=z, cmap=plt.cm.get_cmap('Dark2', 6),s=300)
    axs[i].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    axs[i].set_xlabel(coms[i],fontdict=dict(weight='bold'))
    fig.colorbar(pcm, ax=axs[i],ticks=list(range(int(min(z)),int(max(z)+1),int((max(z)+1)/4))),label = 'gene overlaps')
    
    
plt.rcParams.update({'font.size': 20})
fig.supxlabel(' Adjusted P-value',fontdict=dict(weight='bold')); fig.supylabel('GO Pathways',fontdict=dict(weight='bold'))
fig.suptitle(fname.split('_')[0],fontdict=dict(weight='bold'))

#%%

'''PLOT 3.'''

import matplotlib.ticker as mticker

'''
sort and plot the adjusted p-value distribution among the two partition enrichments, to compare them
'''
cons_pdis = []
coex_pdis = []

for k,v in com_paths.items():
    mu = v[0].head(10)['p-adj'].mean(); sd = v[0].head(10)['p-adj'].std()
    cons_pdis.append([round(mu,4),round(sd,4)])
cons_pdis = list(zip(*cons_pdis)); cons_pdis = pd.DataFrame({'mu':cons_pdis[0],'sd':cons_pdis[1]})
cons_pdis = cons_pdis.sort_values(by=['mu'],ascending=True).reset_index(drop=True)

for k,v in com_paths_coex.items():
    mu = v[0].head(10)['p-adj'].mean(); sd = v[0].head(10)['p-adj'].std()
    coex_pdis.append([round(mu,4),round(sd,4)])
coex_pdis = list(zip(*coex_pdis)); coex_pdis = pd.DataFrame({'mu':coex_pdis[0],'sd':coex_pdis[1]})
coex_pdis = coex_pdis.sort_values(by=['mu'],ascending=True).reset_index(drop=True)


x1 = [z+1 for z in list(range(len(cons_pdis)))] ; y1 = cons_pdis['mu'].tolist(); e1 = cons_pdis['sd'].tolist()
x2 = [z+1 for z in list(range(len(coex_pdis)))] ; y2 = coex_pdis['mu'].tolist(); e2 = coex_pdis['sd'].tolist()

plt.rcParams.update({'font.size': 22})

plt.ticklabel_format(style='plain', axis='x', useOffset=False)

plt.errorbar(x1,y1,yerr=e1,fmt='o',label='consensus'); plt.errorbar(x2,y2,yerr=e2,fmt='o',label='coexpression')
lab = fname.split('_')[0]
plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(1))
plt.gca().set_xlim(left=0)
plt.xticks(rotation=90)
plt.xlabel(lab+' community index'); plt.ylabel('Adjusted P-value')
plt.legend(loc="upper left")

plt.grid(True)
