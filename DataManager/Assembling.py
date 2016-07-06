# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:30:31 2016

@author: huan
"""
import sqlite3
import pandas as pd

def GetMaterialBasic(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of material basic mechanical
    """
    df=pd.read_sql('SELECT * FROM Material_Properties_Basic_Mechanical',conn)
    return df
    
def GetSteel(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of steel
    """
    df=pd.read_sql('SELECT * FROM Material_Properties_Steel',conn)
    return df
    
def GetBeamSections(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of beam sections
    """
    df=GetBeamSections(conn)
    df=pd.read_sql('SELECT * FROM Beam_Prop_General',conn)
    
    return df

def GetNodesCoordiniate(conn):
    df=pd.read_sql('SELECT * FROM Node_Coordinates',conn)
    return df
    
def GetBeams(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of beam
    """
    df=pd.read_sql('SELECT * FROM Connectivity_Beam',conn)
    return df
    
def GetBeamSections(conn):
    """
    conn: sqlite database connection.\n
    return: dataframe of beam sections
    """
    mdf=GetMaterialBasic(conn)
    df1=pd.read_sql('SELECT * FROM Beam_Section_Assignments',conn)
    
    #return EA,EI,GJ,etc.
    
    return df

if __name__=='__main__':
    conn=sqlite3.connect('C:\\huan\\Test\\Test.sqlite')
    try:
        a=GetBeams(conn)
        print(a)
    finally:
        conn.close()
