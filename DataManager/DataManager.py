# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:26:46 2016

@author: huan
"""
import pandas as pd

class ModelData:
    """
    ModelData is used to store data in the RAM with DataFrame
    """
    def __init__(self):
        self.dataFrames={}
        
    def Open(self,path):
        return False
    
    def Save(self,path):
        for df in self.dataFrames:
            df.to_sql()
        

        

    