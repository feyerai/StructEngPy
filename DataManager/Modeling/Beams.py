# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:32:51 2016

@author: huan
"""
import pandas as pd
        
def CreateTable(md):
    """
    md: ModelData
    """
    if 'ConnectivityBeam' not in md.dataFrames.keys():
        md.dataFrames['ConnectivityBeam']=pd.DataFrame({
        'NodeI':[],
        'NodeJ':[]
        })
    if 'FrameSectionAssignments' not in md.dataFrames.keys():
        md.dataFrames['FrameSectionAssignments']=pd.DataFrame({
        'Section':[],
        'MatProp':[]
        })
    
def AddBeams(md,beams):
    """
    md: ModelData\n
    beams: a list of tuples in (name,nodei,nodej,section)
    """
    for beam in beams:
        md.dataFrames['ConnectivityBeam']=md.dataFrames['ConnectivityBeam'].append(pd.Series({
        'NodeI':str(beam[1]),
        'NodeJ':str(beam[2])
        },index=str(beam[0])))
        md.dataFrames['FrameSectionAssignments'].append(pd.Series({
        'Section':str(beam[3]),
        'MatProp':'Default'
        },index=str(beam[0])))

def Count(md):
    """
    md: ModelData
    """
    return md.dataFrames['ConnectivityBeam'].shape[0]

if __name__=='__main__':
    df=pd.read_csv('d:\\testnode.csv',header=None)
    nodes=[]
    for i in range(df.shape[0]):
        node=df.ix[i]
        nodes.append(Node(node[0],node[0],node[1],node[2],node[3],node[4],node[5]))
    dm=Nodes('testsql.sqlite')
    try:
        dm.Connect()
        dm.CreateTable()
        dm.AddCartesian(nodes)
    finally:
        dm.Close()