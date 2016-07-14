# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 09:45:00 2016
Wind load calculation with GB50007-2010
define CalMuzBetaz
define cors=A #Coarse type
Coded by: Zhuoju Huang
E-mail: zhuoju36@hotmail.com
"""
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

def WindLoad(w0,cors,mu_s,f1,h,Ht,Bt,Bb,myB_z,myfi1_z):  
    #####################################################################
    #  Input:                                                           #
    #                   #w0         #Basic wind pressure                #
    #                   #mu_s       #Shape coefficient                  #
    #                   #f1         #1st order frequency                #
    #                   #h          #Height of the measure point        #
    #                   #Ht         #Top height                         #
    #                   #Bt         #Top width                          #
    #                   #Bb         #Bottom width                       #
    #                   #myB_z      #Width at height z                  #
    #                   #myfi1_z    #fi1 at height z                    #
    #  Return:                                                          #
    #                   #w          #windpressure                       #
    #####################################################################

    #Tables======================================================================================================================================================================
    #Table 4.2.6 in GB50009-2010-------------------------------------------------------------------------------------------------------------------------------------------------
    heig=[5,   10,  15,  20,  30,  40,  50,  60,  70,  80,  90,  100, 150, 200, 250, 300, 350, 400, 450, 500, 550]
    mu={
    'A':[1.09,   1.28,   1.42,   1.52,   1.67,   1.79,   1.89,   1.97,   2.05,   2.12,   2.18,   2.23,   2.46,   2.64,   2.78,   2.91,   2.91,   2.91,   2.91,   2.91,   2.91],
    'B':[1.00,   1.00,   1.13,   1.23,   1.39,   1.52,   1.62,   1.71,   1.79,   1.87,   1.93,   2,      2.25,   2.46,   2.63,   2.77,   2.91,   2.91,   2.91,   2.91,   2.91],
    'C':[0.65,   0.65,   0.65,   0.74,   0.88,   1,      1.1,    1.2,    1.28,   1.36,   1.43,   1.5,    1.79,   2.03,   2.24,   2.43,   2.6,    2.76,   2.91,   2.91,   2.91],
    'D':[0.51,   0.51,   0.51,   0.51,   0.51,   0.6,    0.69,   0.77,   0.84,   0.91,   0.98,   1.04,   1.79,   2.03,   2.24,   2.43,   2.6,    2.76,   2.91,   2.91,   2.91]
    }
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #Table 8.4.5-1-----------------------
    k1 =[0.944,   0.67,    0.295,   0.112]
    a11=[0.155,   0.187,   0.261,   0.346]
    k2 =[1.276,   0.910,   0.404,   0.155]
    a12=[0.186,   0.218,   0.292,   0.376]
    #-------------------------------------
    
    #Table 8.4.5-2--------------------------------------------------------------------------
    BRat=[1,     0.9,     0.8,     0.7,     0.6,     0.5,     0.4,     0.3,     0.2,     0.1]
    thetav0=[ 1.00,   1.10,    1.20,    1.32,    1.50,    1.75,    2.08,    2.53,    3.30,    5.60]
    #----------------------------------------------------------------------------------------
    #End Tables======================================================================================================================================================================
    
    gg=2.5
    I10=0.12
    
    xi1=0.01 #steel structure
    
    #Muz=================================================================================================
#    for ii in range(21):
#        if ii==0:
#            if h<heig[ii]:
#                mu_z=mu[cors][ii]
#            elif ii==18:
#                if h>=heig[ii]:
#                    mu_z=mu[cors][ii]
#            else:
#                if h>=heig[ii]:
#                    if h<heig[ii+1]:
#                        mu_z=mu[cors][ii]+mu[cors][ii+1]-mu[cors][ii]*(h-heig[ii])/(heig(ii+1)-heig[ii])
    mu_z=np.interp(h,heig,mu[cors])
    #end muz==============================================================================================
    
    #Betaz================================================================================================
    #Cal kw in 8.4.4-2
    if cors == 'A':
        kw=1.28
    elif cors=='B':
        kw=1.00
    elif cors=='C':
        kw=0.54
    elif cors=='D':
        kw=0.26
        
    x1=30*f1/np.sqrt(kw*w0)
    R=np.sqrt(np.pi*x1**2/6/xi1/(1+x1**2)**(4/3))
    #define k, a1
    if cors=='A':
        k=k2[0]
        a1=a12[0]
    elif cors=='B':
        k=k2[1]
        a1=a12[1]
    elif cors=='C':
        k=k2[2]
        a1=a12[2]
    elif cors=='D':
        k=k[3]
        a1=a12[3]
    rhoz=10*np.sqrt(Ht+60*np.exp(-Ht/60)-60)/Ht
    rhox=10*np.sqrt(Bt+50*np.exp(-Bt/50)-50)/Bt #Not clear if B should be constant
    BGF=k*Ht**a1*rhox*rhoz/mu_z*myfi1_z
    thetav=1 #coordinate factors
    thetab=myB_z/Bb
    beta_z=1+2*gg*I10*BGF*np.sqrt(1+R**2)*thetav*thetab
    #end betaz================================================================= 
    return mu_z*beta_z*w0*mu_s

#nodeFile='C:\\huan\\Project\\1603-Zhaode_T3\\repository\\T3_overall\\windload\\node.txt'
#quadFile='C:\\huan\\Project\\1603-Zhaode_T3\\repository\\T3_overall\\windload\\quad.txt'
#nodeDF=pd.read_table(nodeFile,index_col=0)
#quadDF=pd.read_table(quadFile,index_col=0)
#
#q=pd.DataFrame(columns=['NR','X','Y','Z'])
#for i in quadDF.index:
#    n1=int(quadDF.ix[i]['node1'])
#    n2=int(quadDF.ix[i]['node2'])
#    n3=int(quadDF.ix[i]['node3'])
#    n4=int(quadDF.ix[i]['node4'])
#    x=nodeDF.ix[n1]['X [m]']+nodeDF.ix[n2]['X [m]']+nodeDF.ix[n3]['X [m]']+nodeDF.ix[n4]['X [m]']
#    y=nodeDF.ix[n1]['Y [m]']+nodeDF.ix[n2]['Y [m]']+nodeDF.ix[n3]['Y [m]']+nodeDF.ix[n4]['Y [m]']
#    z=nodeDF.ix[n1]['Z [m]']+nodeDF.ix[n2]['Z [m]']+nodeDF.ix[n3]['Z [m]']+nodeDF.ix[n4]['Z [m]']
#    q=q.append({'NR':i,'X':x/4,'Y':y/4,'Z':z/4},ignore_index=True)
#print(q)
#
#kmeansModel = KMeans(n_clusters=35, random_state=1)
#a=kmeansModel.fit_predict(np.array([q.iloc[:,3].values]).T)
#q['G']=a

try:
    path='C:\\huan\\Project\\1603-Zhaode_T3\\repository\\T3_overall\\'
    file1=open(path+'testwindY_.dat','w+')
    file2=open(path+'testwindX.dat','w+')
    file3=open(path+'testwindY.dat','w+')
    file4=open(path+'testwindX_.dat','w+')
    H=150
    ztoh_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    fi1_list=[0.02,0.08,0.17,0.27,0.38,0.45,0.67,0.74,0.86,1.00]
    for i in q.index:
        w0=0.55
        Bx=40
        By=34
        f1=0.67
        grp=int(q.ix[i]['NR'])//10000
        z=q.ix[i]['Z']
        ztoh=z/H
        fi1 = np.interp(ztoh,ztoh_list,fi1_list)
        mu_s_list=[0.8,0.7,0.5]
        if grp==101:
            w = WindLoad(w0,'A',mu_s_list[0],f1,z,H,By,By,By,fi1)
            file1.write('quad %d'%q.ix[i]['NR']+' type pyy -%6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,By,By,By,fi1)
            file2.write('quad %d'%q.ix[i]['NR']+' type pyy %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[2],f1,z,H,By,By,By,fi1)
            file3.write('quad %d'%q.ix[i]['NR']+' type pyy %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,By,By,By,fi1)
            file4.write('quad %d'%q.ix[i]['NR']+' type pyy %6.5f\n'%w)
        elif grp==102:
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,Bx,Bx,Bx,fi1)
            file1.write('quad %d'%q.ix[i]['NR']+' type pxx -%6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[0],f1,z,H,Bx,Bx,Bx,fi1)
            file2.write('quad %d'%q.ix[i]['NR']+' type pxx  %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,Bx,Bx,Bx,fi1)
            file3.write('quad %d'%q.ix[i]['NR']+' type pxx -%6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[2],f1,z,H,Bx,Bx,Bx,fi1)
            file4.write('quad %d'%q.ix[i]['NR']+' type pxx -%6.5f\n'%w)
        elif grp==103:
            w = WindLoad(w0,'A',mu_s_list[2],f1,z,H,By,By,By,fi1)
            file1.write('quad %d'%q.ix[i]['NR']+' type pyy -%6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,By,By,By,fi1)
            file2.write('quad %d'%q.ix[i]['NR']+' type pyy -%6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[0],f1,z,H,By,By,By,fi1)
            file3.write('quad %d'%q.ix[i]['NR']+' type pyy  %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,By,By,By,fi1)
            file4.write('quad %d'%q.ix[i]['NR']+' type pyy -%6.5f\n'%w)
        elif grp==104:
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,Bx,Bx,Bx,fi1)
            file1.write('quad %d'%q.ix[i]['NR']+' type pxx %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[2],f1,z,H,Bx,Bx,Bx,fi1)
            file2.write('quad %d'%q.ix[i]['NR']+' type pxx %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[1],f1,z,H,Bx,Bx,Bx,fi1)
            file3.write('quad %d'%q.ix[i]['NR']+' type pxx %6.5f\n'%w)
            w = WindLoad(w0,'A',mu_s_list[0],f1,z,H,Bx,Bx,Bx,fi1)
            file4.write('quad %d'%q.ix[i]['NR']+' type pxx -%6.5f\n'%w)
finally:
    file1.close()
    file2.close()
    file3.close()
    file4.close()
    
    