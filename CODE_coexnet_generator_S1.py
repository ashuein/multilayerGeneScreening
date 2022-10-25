#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 05:34:39 2022

@author: giriraj
code to run RloW_PC from python script
"""


import rpy2.robjects as robjects
import pandas as pd
import os
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

robjects.r('''
           
           library(RLowPC)
           library(corpcor)
           library("igraph")
    
    processData <- function(df,fname,writeNet = T)
    {
      path = '/home/giriraj/Documents/ncode/python/paper2_codes/datasets/'
      
     
      rownames(df) <- df$Symbol
      
      #selecting only significant ones
      df <- df[df$FoldChange>=1 | df$FoldChange<=-1,]
      #df <- df[df$padj<=0.1 & df$pvalue<=0.05,] #This filter already applied at the data file generation; request the same from user
      
      df <- subset(df,select = -c(Symbol,FoldChange))
      df <- t(df)
      
      data.exp <- df
      
      genes <- colnames(data.exp)
  
      ##infer PC network
      inf.pcor<-abs(pcor.shrink(data.exp)[1:length(genes),1:length(genes)])
      diag(inf.pcor)<-0
      
      
      ##inf RLowPC
      reduction.sapce<-na.omit(adjmatrix2edgelist(adjmatrix = inf.pcor,directed = F,order = T))
      
      thr <- mean(reduction.sapce$weight) + 1 * sd(reduction.sapce$weight) # thr = mean + 2 * sd
      reduction.sapce <- reduction.sapce[reduction.sapce$weight>=thr,]
      
      inf.RLowPC.edge<-RLowPC(data.exp = data.exp,edgelist = reduction.sapce,
                              method = 'pearson',pc.estimator = 'shrink')
      inf.RLowPC.edge$cor.weight<-abs(inf.RLowPC.edge$cor.weight)
      inf.RLowPC<-edgelist2adjmatrix(inf.RLowPC.edge[,c(1,2,4)],genes = genes,directed = T)
      inf.RLowPC<-abs(inf.RLowPC)
      inf.RLowPC<-pmax(inf.RLowPC,t(inf.RLowPC))
      
      
      
      ##infer first order PC based on reduction space
      ###first PC
      inf.firstPC.edge<-firstPC(data.exp = data.exp,edgelist = reduction.sapce,
                                method = 'pearson',controlist = NULL)
      inf.firstPC<-edgelist2adjmatrix(inf.firstPC.edge[,1:3],genes = genes,directed = F)
      
      
      ###second PC
      #inf.secondPC.edge<-secondPC(data.exp = data.exp,edgelist = reduction.sapce,method = 'pearson',controlist = NULL)
      #inf.secondPC<-edgelist2adjmatrix(inf.secondPC.edge[,1:3],genes = genes,directed = F)
      
      if (writeNet==T) {
        g <- graph.adjacency(inf.firstPC,mode="undirected",weighted=TRUE,diag=FALSE,add.rownames=TRUE)
    
        g <- delete_edges(g, E(g)[which(E(g)$weight==0)])
        g <- delete.vertices(g, degree(g)==0)
        
        fname <- paste(path,fname,sep='')
        write_graph(g,file = fname, format = c("ncol"))
      }
      
      return(inf.firstPC)
    }
    
           ''')
           
           
fnCoexNet = robjects.r['processData']

#%%


os.chdir('~/paper2_codes/') # set the current directory here
cwd = os.getcwd()

read_csv = robjects.r['read.csv']


#NOTE: To process individual files: un-comment the file to be processed and comment out the rest of the names.

#fname = 'modified_sod_expData.csv'; fname1 = 'sod_coex_net.csv'
#fname = 'modified_tdp43_expData.csv'; fname1 = 'tdp43_coex_net.csv'
#fname = 'GSE145677_DEGs.csv'; fname1 = 'GSE145677_coex_net.csv'
#fname = 'GSE158264_DEGs.csv'; fname1 = 'GSE158264_coex_net.csv'
fname = 'GSE40438_DEGs.csv'; fname1 = 'GSE40438_coex_net.csv'



def coex_net_generator(fname,cwd,coexfile):
    df = pd.read_csv(cwd+'/datasets/'+fname).dropna(subset=['Symbol'])
    df['Symbol'] = df['Symbol'].apply(lambda x:x.split(';')[0])
    df = df.dropna().reset_index(drop=True)
    df = df.groupby(by='Symbol',as_index=False).mean()
    if 'ID' in df.columns: df = df.drop('ID', 1)
    
    try:
        with localconverter(ro.default_converter + pandas2ri.converter): #conversion Pandas to R dataframe
            df = ro.conversion.py2rpy(df)
            
        val = fnCoexNet(df,coexfile)
        return val
    except:
        df = df[~df['Symbol'].str.contains("LOC")]
        with localconverter(ro.default_converter + pandas2ri.converter): #conversion Pandas to R dataframe
            df = ro.conversion.py2rpy(df)
            
        val = fnCoexNet(df,coexfile)
        return val
        
val = coex_net_generator(fname,cwd,fname1)




















           
           
