# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 22:13:07 2016

@author: HZJ
"""
import DataManager.Definition.Properties.SectionParameters as sp
import pandas as pd

def CreateTable(md):
    """
    md: ModelData
    """
    if 'BeamPropGeneral' not in md.dataFrames.keys():
        md.dataFrames['BeamPropGeneral']=pd.DataFrame({
        'Name':[],
        'Material':[],
        'Shape':[],
        't3':[],
        't2':[],
        'tf':[],
        'tw':[],
        'Area':[],
        'TorsConst':[],
        'I33':[],
        'I22':[],
        'AS2':[],
        'AS3':[],
        'S33':[],
        'S22':[],
        'Z33':[],
        'Z22':[],
        'R33':[],
        'R22':[],
        'Color':[],
        'FromFile':[],
        'AMod':[],
        'A2Mod':[],
        'A3Mod':[],
        'JMod':[],
        'I2Mod':[],
        'I3Mod':[],
        'MMod':[],
        'WMod':[],
        'GUID':[],
        'Notes':[]
        })

def AddQuick(md,material,profile,commit=True):
    """
    md: ModelData
    profile: format like H400x200x10x12
    now only I-profile is allowable
    """
    if profile[0]=='H':
        size=profile[1:].split('x')
        h=eval(size[0])
        b=eval(size[1])
        tw=eval(size[2])
        tf=eval(size[3])
        section=sp.IProfile(h,b,tw,tf)
    
    md.dataFrames['BeamPropGeneral'].append({
        'Name':profile,
        'Material':material,
        'Shape':'I',
        't3':size[0],
        't2':size[1],
        'tf':size[3],
        'tw':size[2],
        'Area':section.A,
        'TorsConst':section.J,
        'I33':section.I33,
        'I22':section.I22,
        'AS2':0,
        'AS3':0,
        'S33':0,
        'S22':0,
        'Z33':0,
        'Z22':0,
        'R33':0,
        'R22':0,
        'Color':0x000000,
        'FromFile':0,
        'AMod':0,
        'A2Mod':0,
        'A3Mod':0,
        'JMod':0,
        'I2Mod':0,
        'I3Mod':0,
        'MMod':0,
        'WMod':0,
        'GUID':'0',
        'Notes':0
        })
    
