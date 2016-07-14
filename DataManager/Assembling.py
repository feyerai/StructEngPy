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
    
def GetSteel(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of steel
    """
    df=md.dataFrames['MaterialPropertiesSteel']
    return df
    
def GetBeamSections(conn):
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
    
def GetBeams(conn):
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
    mdf=GetMaterialBasic(md)
    bsdf=GetBeamSections(md)
    EA=[]
    EI33=[]
    EI22=[]
    GJ=[]
    for i in range(bsdf.shape[0]):
        E=mdf[mdf.Name==sdf.ix[i]['Material']]['E1'].value
        G=mdf[mdf.Name==sdf.ix[i]['Material']]['G12'].value
        A  =sdf.ix[i]['Area']
        I33=sdf.ix[i]['I33']
        I22=sdf.ix[i]['I22']
        J  =sdf.ix[i]['TorsConst']
        EA.append(E*A)
        EI33.append(E*I33)
        EI22.append(E*I22)
        GJ.append(G*J)
    
    sdf=pd.DataFrame()
    sdf['EA']=pd.Series(EA,index=bsdf.index)
    sdf['EI33']=pd.Series(EI33,index=bsdf.index)    
    sdf['EI22']=pd.Series(EI22,index=bsdf.index)
    sdf['GJ']=pd.Series(GJ,index=bsdf.index)
    
    cdf=md.dataFrames['ConnectivityBeam']
    bsdf=md.dataFrames['BeamSectionAssignments']
    
    ndf=md.dataFrames['NodeCoordinates']
    for i in range(cdf.shape[0]):
        n1=ndf[ndf.Name==sdf.ix[i]]
        n2=ndf[ndf.Name==sdf.ix[i]]
        length.append(np.sqrt((n1['X'].value-n2['X'].value)**2+(n1['Y'].value-n2['Y'].value)**2+(n1['Y'].value-n2['Y'].value)**2))
    
    df=pd.DataFrame(cdf['Name'])
    df['Length']=pd.Series(length,index=cdf.index)
    for i in range(bsdf.shape[0]):
        
    
    
    
    #return EA,EI,GJ,etc.
    
    return df

if __name__=='__main__':
    conn=sqlite3.connect('C:\\huan\\Test\\Test.sqlite')
    try:
        a=GetBeams(conn)
        print(a)
    finally:
        conn.close()
