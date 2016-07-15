# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:26:46 2016

@author: huan
"""
import shutil
import sqlite3
import pandas as pd

class ModelData:
    """
    ModelData is used to store data in the RAM with DataFrame
    """
    def __init__(self):
        self.dataFrames={}
        
    def Open(self,file):
        return False
    
    def Save(self,db):
        tempPath='F:\\test\\'
        temp=tempPath+'temp.sqlite'
        try:
            conn=sqlite3.connect(temp)
            for name in self.dataFrames.keys():
                self.dataFrames[name].to_sql(name,conn)              
        finally:
            conn.close()
        shutil.move(temp,db)