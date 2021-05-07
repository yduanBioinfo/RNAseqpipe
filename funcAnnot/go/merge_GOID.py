#!/usr/bin/env python3

""" Add GOID column from GOMFID, GOBPID and GOCCID
    GOBPID  GOCCID  GOMFID  Pvalue  OddsRatio       ExpCount        Count   Size    Term    GODomain        Direction       left    right   event_type
    GO:0000209                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein polyubiquitination      BP      over    T4    T5      All
    GO:0070936                      0.009489876710779171    418.0   0.00950118764845606     1       2       protein K48-linked ubiquitination       BP      over    T4      T5      All
    GO:0006414                      0.0142178486596539      208.5   0.0142517814726841      1       3       translational elongation        BP      over    T4    T5      All
        GO:0005765              0.0169488797321661      174.5   0.0169971671388102      1       3       lysosomal membrane      CC      over    T4      T5      All
"""
import sys
import pandas as pd

df = pd.read_csv(sys.stdin,sep="\t",header=0)
df['GOID']='N'
df.loc[df['GODomain']=='BP','GOID'] = df.loc[df['GODomain']=='BP','GOBPID']
df.loc[df['GODomain']=='MF','GOID'] = df.loc[df['GODomain']=='MF','GOMFID']
df.loc[df['GODomain']=='CC','GOID'] = df.loc[df['GODomain']=='CC','GOCCID']
df.to_csv(sys.stdout,sep="\t",index=False)
