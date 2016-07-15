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
    if 'NodeCoordinates' not in md.dataFrames.keys():
        md.dataFrames['NodeCoordinates']=pd.DataFrame({
        'CoordSys':[],
        'CoordType':[],
        'X':[],
        'Y':[],
        'Z':[]
        })
    
def AddCartesian(md,nodes):
    """
    md: ModelData
    nodes: a list of tuples in (name,x,y,z)
    """
    for node in nodes:
        md.dataFrames['NodeCoordinates']=md.dataFrames['NodeCoordinates'].append(pd.DataFrame({
        'CoordSys':'Global',
        'CoordType':'Cartesian',
        'X':node[1],
        'Y':node[2],
        'Z':node[3]
        },index=[str(node[0])]))
            
def Count(md):
    """
    md: ModelData
    """
    return md.dataFrames['NodeCoordinates'].shape[0]
    
def GetCoordCartesian():
    return False        

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

