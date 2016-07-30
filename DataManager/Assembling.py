# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:30:31 2016

@author: huan

This package preprocess the relative parameters for analysis.
"""
import os,sys 
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) 
sys.path.insert(0,parentdir)

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
    df=db.dataFrames['MaterialPropertiesSteel']
    return df
    
def GetBeamSectionProp(db):
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


    

    
def GetBeamStiffness(db):
    """
    md: ModelData object..\n
    return: dataframe of beam stifness info with HID
    """
    mdf=GetMaterialBasic(db) #material
    bsdf=GetBeamSectionProp(db) #section
    badf=GetBeamSectionAssignments(db)#section assignment
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
        bEA.append(sdf['EA'][badf['Section'][idx]])
        bEI33.append(sdf['EI33'][badf['Section'][idx]])
        bEI22.append(sdf['EI22'][badf['Section'][idx]])
        bGJ.append(sdf['GJ'][badf['Section'][idx]])
    df=cdf.copy()
    df['L']=pd.Series(length,index=cdf.index)
    df['EA']=pd.Series(bEA,index=cdf.index)
    df['EI33']=pd.Series(bEI33,index=cdf.index)
    df['EI22']=pd.Series(bEI22,index=cdf.index)
    df['GJ']=pd.Series(bGJ,index=cdf.index)
    return df

if __name__=='__main__':
    node=GetBeams('F:\\Test\\Test.sqlite')
    print(node)
    