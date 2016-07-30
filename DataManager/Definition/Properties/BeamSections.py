# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 22:13:07 2016

@author: HZJ
"""
import numpy as np
import pandas as pd

def CreateTable(md):
    """
    md: ModelData
    """
    if 'BeamPropGeneral' not in md.dataFrames.keys():
        md.dataFrames['BeamPropGeneral']=pd.DataFrame({
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
        section=IProfile(h,b,tw,tf)
    CreateTable(md)
    md.dataFrames['BeamPropGeneral']=md.dataFrames['BeamPropGeneral'].append(pd.DataFrame({
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
        },index=[profile]))


class Section:
    def __init__(self,A,I33,I22,W33,W22):
        """
        A   - [mm]\n
        I33 - [mm]\n
        I22 - [mm]\n
        W33 - [mm]\n
        W22 - [mm]\n
        """
        self.A=A
        self.I33=I33
        self.I22=I22
        self.W33=W33
        self.W22=W22
        self.i33=np.sqrt(self.I33/self.A)
        self.i22=np.sqrt(self.I22/self.A)
        
class Rectangle(Section):
    def __init__(self,h,b):
        """
        h - [mm]\n
        b - [mm]\n
        """
        self.h=h
        self.b=b
        A=h*b
        I33=b*h**3/12
        I22=h*b**3/12
        W33=I33/h*2
        W22=I22/b*2
        Section.__init__(self,A,I33,I22,W33,W22)
#        self.gamma33=1.05
#        self.gamma22=1.05


class Circle(Section):
    def __init__(self,d):
        """
        d - [mm]
        """
        self.d=d
        A=np.pi*d**2/4
        I33=np.pi*d**4/64
        I22=I33
        W33=I33/d*2
        W22=W33
        Section.__init__(self,A,I33,I22,W33,W22)
#        self.gamma33=1.15
#        self.gamma22=1.15
        
class Tube(Section):
    def __init__(self,d,t,fab='r'):
        """
        d - [mm]\n
        t - [mm]\n
        fab - fabrication\n
            'r' - rolled\n
            'w' - welded\n
        """
        self.d=d
        self.t=t
        A=np.pi*d**2/4-np.pi*(d-2*t)**2/4
        I33=np.pi*d**4/64*(1-((d-2*t)/d)**4)
        I22=I33
        W33=I33/d*2
        W22=W33
        Section.__init__(self,A,I33,I22,W33,W22)
        self.gamma33=1.15
        self.gamma22=1.15
        if fab=='r':
            self.cls33='b'
            self.cls22='b'
        elif fab=='w':
            self.cls33='c'
            self.cls22='c'
        else:
            raise ValueError('wrong fabrication!')

class HollowBox(Section):
    def __init__(self,h,b,tw,tf,fab='r'):
        """
        h - [mm]\n
        b - [mm]\n
        tw- [mm]\n
        tf- [mm]\n
        fab - fabrication\n
            'r' - rolled\n
            'w' - welded\n
        """
        self.h=h
        self.b=b
        self.tw=tw
        self.tf=tf
        A=h*b-(h-2*tf)*(b-2*tw)
        I33=b*h**3/12-(b-2*tw)*(h-2*tf)**3/12
        I22=h*b**3/12-(h-2*tf)*(b-2*tw)**3/12
        W33=I33/h*2
        W22=I22/b*2
        Section.__init__(self,A,I33,I22,W33,W22)
        self.gamma33=1.05
        self.gamma22=1.05
        self.cls33='c'
        self.cls22='c'
        
class IProfile(Section):
    def __init__(self,h,b,tw,tf,fab='r'):
        """
        h - [mm]\n
        b - [mm]\n
        tw- [mm]\n
        tf- [mm]\n
        fab - fabrication\n
            'r' - rolled\n
            'w' - welded\n
        """
        self.h=h
        self.b=b
        self.tw=tw
        self.tf=tf
        A=b*tf*2+tw*(h-2*tf)
        I33=b*h**3/12-(b-tw)*(h-2*tf)**3/12
        I22=2*tf*b**3/12+(h-2*tf)*tw**3/12
        W33=I33/h*2
        W22=I22/b*2
        self.J=1/3*(b*tf**4*2+(h-2*tf)*tw**4)
        Section.__init__(self,A,I33,I22,W33,W22)
        self.gamma33=1.05
        self.gamma22=1.2
        self.cls33='c'
        self.cls22='c'