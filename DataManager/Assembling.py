# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:30:31 2016

@author: huan

This package preprocess the relative parameters for analysis.
"""
import sqlite3
import pandas as pd
import numpy as np

def GetMaterialBasic(md):
    """
    md: ModelData object.\n
    return: dataframe of material basic mechanical
    """
    df=md.dataFrames['MaterialPropertiesBasicMechanical']
    return df
    
def GetSteel(md):
    """
    conn: sqlite database connection.\n
    return: dataframe of steel
    """
    df=md.dataFrames['MaterialPropertiesSteel']
    return df
    
def GetBeamSections(md):
    """
    conn: sqlite database connection.\n
    return: dataframe of beam sections
    """
    df=md.dataFrames['BeamPropGeneral']    
    return df

def GetNodesCoordiniate(md):
    """
    md: ModelData object.
    """
    df=md.dataFrames['NodeCoordinates']
    return df
    
def GetBeams(md):
    """
    md: ModelData object.\n
    return: dataframe of beam
    """
    df=md.dataFrames['ConnectivityBeam']
    return df
    
def GetBeamStiffness(md):
    """
    md: ModelData object..\n
    return: dataframe of beam sections
    """
    mdf=GetMaterialBasic(md) #material
    bsdf=GetBeamSections(md) #section
    EA=[]
    EI33=[]
    EI22=[]
    GJ=[]
    for idx,row in bsdf.iterrows():
        E=mdf.ix[row['Material']]['E1'].get_value(0)
        G=mdf.ix[row['Material']]['G12'].get_value(0)
        A  =bsdf.ix[idx]['Area']
        I33=bsdf.ix[idx]['I33']
        I22=bsdf.ix[idx]['I22']
        J  =bsdf.ix[idx]['TorsConst']
        EA.append(E*A)
        EI33.append(E*I33)
        EI22.append(E*I22)
        GJ.append(G*J)
    
    sdf=pd.DataFrame() #section stiffness
    sdf['EA']=pd.Series(EA,index=bsdf.index)
    sdf['EI33']=pd.Series(EI33,index=bsdf.index)    
    sdf['EI22']=pd.Series(EI22,index=bsdf.index)
    sdf['GJ']=pd.Series(GJ,index=bsdf.index)
    
    cdf=md.dataFrames['ConnectivityBeam']
    ndf=md.dataFrames['NodeCoordinates']
    length=[]
    for idx,row in cdf.iterrows():
        n1=ndf.ix[row['Joint1']]
        n2=ndf.ix[row['Joint2']]
        length.append(np.sqrt((n1['X'].value-n2['X'].value)**2+(n1['Y'].value-n2['Y'].value)**2+(n1['Y'].value-n2['Y'].value)**2))
    
    df=pd.DataFrame({'L':length},index=cdf.index)
    df['EA']=sdf['EA']
    df['EI33']=sdf['EI33']
    df['EI22']=sdf['EI22']
    df['GJ']=sdf['GJ'] 
    return df
    