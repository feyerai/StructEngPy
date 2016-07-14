# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:32:51 2016

@author: huan
"""
import pandas as pd

#def CreateTable(conn,commit=True):
#    sqls=['CREATE TABLE IF NOT EXISTS Connectivity_Beam(\
#    Name TEXT UNIQUE PRIMARY KEY, \
#    NodeI TEXT, \
#    NodeJ TEXT\
#    )']
#    sqls.append('CREATE TABLE IF NOT EXISTS Frame_Section_Assignments(\
#    Name TEXT UNIQUE PRIMARY KEY, \
#    Section TEXT, \
#    MatProp TEXT\
#    )')
#    cu=conn.cursor()
#    for sql in sqls:
#        cu.execute(sql)
#    cu.close
#    if commit:
#        conn.commit()
#    
#def AddBeams(conn,beams,commit=True):
#    """
#    conn: sqlite database connection\n
#    beams: a list of tuples in (name,nodei,nodej,section)
#    """
#    #Check nodes and sections first!
#    cu=conn.cursor()
#    sql1='INSERT INTO Connectivity_Beam VALUES (?,?,?)'
#    sql2='INSERT INTO Beam_Section_Assignments VALUES (?,?,?)'
#    args1=[]
#    args2=[]
#    for beam in nodes:
#        args1.append((str(beam[0]),str(beam[1]),str(beam[2])))
#        args2.append((str(beam[0]),str(beam[3]),'Default'))
#    cu.executemany(sql1,args1)    
#    cu.executemany(sql2,args2)
#    cu.close
#    if commit:
#        conn.commit()
        
def CreateTable(md):
    """
    md: ModelData
    """
    if 'ConnectivityBeam' not in md.dataFrames.keys():
        md.dataFrames['ConnectivityBeam']=pd.DataFrame({
        'Name':[],
        'NodeI':[],
        'NodeJ':[]
        })
    if 'FrameSectionAssignments' not in md.dataFrames.keys():
        md.dataFrames['FrameSectionAssignments']=pd.DataFrame({
        'Name':[],
        'Section':[],
        'MatProp':[]
        })
    
def AddBeams(md,beams):
    """
    md: ModelData\n
    beams: a list of tuples in (name,nodei,nodej,section)
    """
    for beam in beams:
        md.dataFrames['ConnectivityBeam'].append({
        'Name':str(beam[0]),
        'NodeI':str(beam[1]),
        'NodeJ':str(beam[2])
        })
        md.dataFrames['FrameSectionAssignments'].append({
        'Name':str(beam[0]),
        'Section':str(beam[3]),
        'MatProp':'Default'
        })
                    
def Count(md):
    """
    md: ModelData
    """
    return md.dataFrames['ConnectivityBeam'].shape[0]
    
        
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