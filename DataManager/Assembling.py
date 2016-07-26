# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:30:31 2016

@author: huan

This package preprocess the relative parameters for analysis.
"""
import sqlite3
import pandas as pd
import numpy as np

def GetMaterialBasic(db):
    """
    db: sqlite data base\n
    return: dataframe of material basic mechanical
    """
    try:
        conn=sqlite3.connect(db)
        df=pd.read_sql('SELECT * FROM MaterialPropertiesBasicMechanical',conn,index_col='index')
    except Exception as e:
        print(e)
        return None
    finally:
        conn.close()
    return df
    
def GetSteel(db):
    """
    db: sqlite data base\n
    return: dataframe of steel
    """
    df=md.dataFrames['MaterialPropertiesSteel']
    return df
    
def GetBeamSections(db):
    """
    db: sqlite data base\n
    return: dataframe of beam sections
    """  
    try:
        conn=sqlite3.connect(db)
        df=pd.read_sql('SELECT * FROM BeamPropGeneral',conn,index_col='index')
    except Exception as e:
        print(e)
        return None
    finally:
        conn.close()
    return df

def GetNodes(db):
    """
    db: sqlite data base\n
    return: dataframe of node coordinates with HID
    """
    try:
        conn=sqlite3.connect(db)
        df=pd.read_sql('SELECT * FROM NodeCoordinates',conn,index_col='index')
        df['HID']=range(df.shape[0])
    except Exception as e:
        print(e)
        return None
    finally:
        conn.close()
    return df
    
def GetBeams(db):
    """
    db: sqlite data base\n
    return: dataframe of beam
    """
    try:
        conn=sqlite3.connect(db)
        df=pd.read_sql('SELECT * FROM ConnectivityBeam',conn,index_col='index')
    except Exception as e:
        print(e)
        return None
    finally:
        conn.close()
    return df
    
def GetBeamStiffness(db):
    """
    md: ModelData object..\n
    return: dataframe of beam stifness info with HID
    """
    mdf=GetMaterialBasic(db) #material
    bsdf=GetBeamSections(db) #section
    cdf=GetBeams(db) #beam connectivity
    ndf=GetNodes(db) #node coordinates
    EA=[]
    EI33=[]
    EI22=[]
    GJ=[]
    for idx,row in bsdf.iterrows():
        E=mdf['E1'][row['Material']]
        G=mdf['G12'][row['Material']]
        A  =bsdf['Area'][idx]
        I33=bsdf['I33'][idx]
        I22=bsdf['I22'][idx]
        J  =bsdf['TorsConst'][idx]
        EA.append(E*A)
        EI33.append(E*I33)
        EI22.append(E*I22)
        GJ.append(G*J)
    
    sdf=pd.DataFrame() #section stiffness
    sdf['EA']=pd.Series(EA,index=bsdf.index)
    sdf['EI33']=pd.Series(EI33,index=bsdf.index)    
    sdf['EI22']=pd.Series(EI22,index=bsdf.index)
    sdf['GJ']=pd.Series(GJ,index=bsdf.index)
    
    length=[]
    bEA=[]
    bEI33=[]
    bEI22=[]
    bGJ=[]
    for idx,row in cdf.iterrows():
        n1=ndf.ix[row['NodeI']]
        n2=ndf.ix[row['NodeJ']]
        length.append(np.sqrt((n1['X']-n2['X'])**2+(n1['Y']-n2['Y'])**2+(n1['Y']-n2['Y'])**2))
        bEA.append()
    
    df=cdf.copy()
    df['L']=pd.Series(length,index=cdf.index)
    for sec in bsdf[]:
    df['EA']=sdf['EA']
    df['EI33']=sdf['EI33']
    df['EI22']=sdf['EI22']
    df['GJ']=sdf['GJ']
    df['HID']=range(df.shape[0])
    return df

if __name__=='__main__':
    mdf=GetBeamStiffness('F:\\Test\\Test.sqlite')    
    