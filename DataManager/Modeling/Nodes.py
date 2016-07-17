# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:32:51 2016

@author: huan
"""

import pandas as pd
    
def AddCartesian(md,nodes):
    """
    md: ModelData
    nodes: a list of tuples in (name,x,y,z)
    """
    if 'NodeCoordinates' not in md.dataFrames.keys():
        md.dataFrames['NodeCoordinates']=pd.DataFrame({
        'CoordSys':[],
        'CoordType':[],
        'X':[],
        'Y':[],
        'Z':[]
        })        
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
    
def DeleteConstraint():
    return False
    
def DeleteLoadDisp():
    return False
    
def DeleteLoadForce():
    return False
    
def DeleteMass():
    return False
    
def DeleteRestraint():
    return False

def DeleteSpring():
    return False
    
def SetConstraint():
    return False

def SetLoadDisp():
    return False

def SetLoadForce(md,node,lc,load):
    """
    node: index of node\n
    lc: name of load case\n
    load: list contains force P1,P2,P3,M1,M2,M3.
    """
    if 'NodalLoads_Force' not in md.dataFrames.keys():
        md.dataFrames['NodalLoads_Force']=pd.DataFrame({
        'Node':[],
        'LC':[],
        'P1':[],
        'P2':[],
        'P3':[],
        'M1':[],
        'M2':[],
        'M3':[]
        }) 
    md.dataFrames['NodalLoads_Force']=md.dataFrames['NodalLoads_Force'].append({
        'Node':str(node),
        'LC':str(lc),
        'P1':load[0],
        'P2':load[1],
        'P3':load[2],
        'M1':load[3],
        'M2':load[4],
        'M3':load[5]
        })
    
def SetLocalAxis():
    return False
    
def SetMass():
    return False
    
def SetSpring():
    return False
        
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

