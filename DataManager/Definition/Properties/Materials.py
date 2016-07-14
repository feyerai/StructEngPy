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
    if 'MaterialPropertiesGeneral' not in md.dataFrames.keys():
         md.dataFrames['MaterialPropertiesGeneral']=pd.DataFrame({
         'Name':[],
         'Type':[],
         'SymType':[],
         'TempDepend':[],
         'Color':[]
         })
         
    if 'MaterialPropertiesBasicMechanical' not in md.dataFrames.keys():
        md.dataFrames['MaterialPropertiesBasicMechanical']=pd.DataFrame({
        'Name':[],
        'UnitWeight':[],
        'UnitMass' :[],
        'E1':[],
        'G12':[],
        'U12':[],
        'A1':[]
        })
        
    if 'MaterialPropertiesSteel' not in md.dataFrames.keys():
        md.dataFrames['MaterialPropertiesSteel']=pd.DataFrame({
        'Name':[],
        'Fy':[],
        'Fu':[],
        'EffFy':[],
        'EffFu':[],
        'SSCurveOpt':[], 
        'SSHysType':[],
        'SHard':[],
        'SMax':[],
        'SRup':[],
        'FinalSlope':[],
        })

def AddQuick(md,std,grade):
    """
    md: ModelData
    std： GB as Chinese standard
    """
    if std=='GB':
        if grade=='Q345':
            md.dataFrames['MaterialPropertiesGeneral'].append({
            'Name':'Q345',
            'Type':'Steel',
            'SymType':'Iso',
            'TempDepend':False,
            'Color':0x000000
            })
            md.dataFrames['MaterialPropertiesBasicMechanical'].append({
            'Name':'Q345',
            'UnitWeight':774.9,
            'UnitMass' :7849,
            'E1':2.000E11,
            'G12':0.3,
            'U12':0.3,
            'A1':1.17e-5
            })
            md.dataFrames['MaterialPropertiesSteel'].append({
            'Name':'Q345',
            'Fy':345,
            'Fu':460,
            'EffFy':0,
            'EffFu':0,
            'SSCurveOpt':0, 
            'SSHysType':0,
            'SHard':0,
            'SMax':0,
            'SRup':0,
            'FinalSlope':0,
            })
            
            
