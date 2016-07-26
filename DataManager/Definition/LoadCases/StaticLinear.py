# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 22:54:09 2016

@author: HZJ
"""
import pandas as pd

def SetCase(md,case):
    """
    md: ModelData\n
    cases: name of load case
    """
    loadtype='Static'
    plc='None'
    modalCase=''
    baseCase=''
    if 'LCDefinintions' not in md.dataFrames.keys():
        md.dataFrames['LCDefinintions']=pd.DataFrame({
        'Type':[],
        'InitialCond':[],
        'ModalCase':[],
        'BaseCase':[],
        'DesTypeOpt':[],
        'DesignType':[],
        'AutoType':[],
        'RunCase':[],
        'GUID':[],
        'Notes':[]
        })
        
    md.dataFrames['LCDefinintions']=md.dataFrames['LCDefinintions'].append(pd.DataFrame({
        'Type':loadtype,
        'InitialCond':plc,
        'ModalCase':modalCase,
        'BaseCase':baseCase,
        'DesTypeOpt':'',
        'DesignType':'',
        'AutoType':'',
        'RunCase':'',
        'GUID':'',
        'Notes':''
        },index=[case]))
    
def SetLoads(md,name,loadtype):
    return False

def SetInitialCase(md,name,iniCase):
    return False
    
def GetInitialCase(md):
    return False
    
def GetLoads(md):
    return False