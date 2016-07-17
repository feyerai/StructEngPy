# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:32:51 2016

@author: huan
"""
import pandas as pd
    
def AddBeams(md,beams):
    """
    md: ModelData\n
    beams: a list of tuples in (name,nodei,nodej,section)
    """
    if 'ConnectivityBeam' not in md.dataFrames.keys():
        md.dataFrames['ConnectivityBeam']=pd.DataFrame({
        'NodeI':[],
        'NodeJ':[]
        })
    if 'BeamSectionAssignments' not in md.dataFrames.keys():
        md.dataFrames['BeamSectionAssignments']=pd.DataFrame({
        'Section':[],
        'MatProp':[]
        })
    for beam in beams:
        md.dataFrames['ConnectivityBeam']=md.dataFrames['ConnectivityBeam'].append(pd.DataFrame({
        'NodeI':str(beam[1]),
        'NodeJ':str(beam[2])
        },index=[str(beam[0])]))
        md.dataFrames['BeamSectionAssignments']=md.dataFrames['BeamSectionAssignments'].append(pd.DataFrame({
        'Section':str(beam[3]),
        'MatProp':'Default'
        },index=[str(beam[0])]))

def Count(md):
    """
    md: ModelData
    """
    return md.dataFrames['ConnectivityBeam'].shape[0]
    
def Delete(md,index):
    return False
    
def DeleteLoadForce():
    return False
    
def DeleteLoadStrain():
    return False
    
def DeleteLoadTemperature():
    return False
    
def DeleteMass():
    return False

def SetLoadForce(md,beam,lc,loadI,loadJ):
    """
    beam: index of beam.\n
    lc: load case\n
    loadI,loadJ: load values indicates the distributed load of P1,P2,P3,M1,M2,M3
    """
    if 'BeamLoadForce' not in md.dataFrames.keys():
        md.dataFrames['BeamLoadForce']=pd.DataFrame({
        'LC':[],
        'P1I':[],
        'P2I':[],
        'P3I':[],
        'M1I':[],
        'M2I':[],
        'M3I':[],
        'P1J':[],
        'P2J':[],
        'P3J':[],
        'M1J':[],
        'M2J':[],
        'M3J':[],
        })
    md.dataFrames['BeamLoadForce']=md.dataFrames['BeamLoadForce'].append(pd.DataFrame({
        'LC':str(lc),
        'P1I':loadI[0],
        'P2I':loadI[1],
        'P3I':loadI[2],
        'M1I':loadI[3],
        'M2I':loadI[4],
        'M3I':loadI[5],
        'P1J':loadJ[0],
        'P2J':loadJ[1],
        'P3J':loadJ[2],
        'M1J':loadJ[3],
        'M2J':loadJ[4],
        'M3J':loadJ[5],
        },index=[str(beam)])
    )

def SetLoadStrain():
    return False
    
def SetLoadTemperature():
    return False
    
def SetLocalAxes():
    return False
    
def SetMass():
    return False
    
def SetReleases(md,beam,resI,resJ):
    """
    beam: index of beam.\n
    resI,resJ: Boolean values indicates the end releases of P,V2,V3,T,M2,M3
    """
    if 'BeamReleaseGeneral' not in md.dataFrames.keys():
        md.dataFrames['BeamReleaseGeneral']=pd.DataFrame({
        'PI':[],
        'V2I':[],
        'V3I':[],
        'TI':[],
        'M2I':[],
        'M3I':[],
        'PJ':[],
        'V2J':[],
        'V3J':[],
        'TJ':[],
        'M2J':[],
        'M3J':[],
        })
    md.dataFrames['BeamReleaseGeneral']=md.dataFrames['BeamReleaseGeneral'].append(pd.DataFrame({
        'PI':resI[0],
        'V2I':resI[1],
        'V3I':resI[2],
        'TI':resI[3],
        'M2I':resI[4],
        'M3I':resI[5],
        'PJ':resJ[0],
        'V2J':resJ[1],
        'V3J':resJ[2],
        'TJ':resJ[3],
        'M2J':resJ[4],
        'M3J':resJ[5],
        },index=[str(beam)])
    )
    
def SetSection():
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