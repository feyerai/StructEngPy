# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:32:51 2016

@author: huan
"""

import pandas as pd

def CreateTable(conn,commit=True):
    cu=conn.cursor()
    sqls=[]
    sqls.append(
    'CREATE TABLE IF NOT EXISTS Material_Properties_General (\
    Name TEXT UNIQUE PRIMARY KEY, \
    Type TEXT, \
    SymType TEXT, \
    TempDepend TEXT, \
    Color REAL)'
    )
    sqls.append(
    'CREATE TABLE IF NOT EXISTS Material_Properties_Basic_Mechanical (\
    Name TEXT UNIQUE PRIMARY KEY, \
    UnitWeight REAL, \
    UnitMass REAL, \
    E1 REAL, \
    G12 REAL, \
    U12 REAL, \
    A1 REAL)'
    )
    sqls.append(
    'CREATE TABLE IF NOT EXISTS Material_Properties_Steel (\
    Name TEXT UNIQUE PRIMARY KEY, \
    Fy REAL, \
    Fu REAL, \
    EffFy REAL, \
    EffFu REAL, \
    SSCurveOpt TEXT, \
    SSHysType TEXT, \
    SHard REAL, \
    SMax REAL, \
    SRup REAL, \
    FinalSlope REAL)'
    )
    for sql in sqls:
        cu.execute(sql)
    cu.close
    if commit:
        conn.commit()

def AddQuick(conn,std,grade,commit=True):
    """
    std： GB as Chinese standard
    """
    cu=conn.cursor()
    if std=='GB':
        if grade=='Q345':
            mpg=('Q345','Steel','Iso','No','0x000000')
            mpb=('Q345','774.9','7849','2.000E11','0.3','0.3','1.17e-5')
            mps=('Q345','345','460','0','0','0','0','0','0','0','0')
            cu.execute('INSERT INTO Material_Properties_General VALUES (?,?,?,?,?)',mpg)
            cu.execute('INSERT INTO Material_Properties_Basic_Mechanical VALUES (?,?,?,?,?,?,?)',mpb)
            cu.execute('INSERT INTO Material_Properties_Steel VALUES (?,?,?,?,?,?,?,?,?,?,?)',mps)
    cu.close
    if commit:
        conn.commit()
        
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
            
            
