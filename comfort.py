# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 17:53:11 2016

@author: huan
"""
import math
w_=6.3
L=6.36*2
C=1
B=C*L
fn=4.72
p0=0.3
beta=0.035
w=w_*B*L
g=9.81
Fp=p0*math.exp(-0.35*fn)
ap=Fp/beta/w*g
beta=0.035

print('C = %.1f'%C)
print('L = %.3f'%L)
print('B = %.3f (A.0.3-2)'%B)
print('w_ = %.3f'%w_)
print('w = %.3f (A.0.3-1)'%w)
print('beta = %.3f'%beta)
print('p0 = %.3f (A.0.2)'%p0)
print('Fp = %.3f (A.0.2-2)'%Fp)
print('ap = %.3f (A.0.2-1)'%ap)