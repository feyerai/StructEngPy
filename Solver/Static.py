# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 14:11:14 2016

@author: HZJ
"""
import numpy as np
import scipy.sparse.linalg as sl

def SolveLinear(aModel,lc,path):
    """
    linear solve the analysis model
    aModel: analysis model
    lc: load case to be solved
    path: path of the matrix files
    """
    if not aModel.isAssembled():
        raise Exception('Not assemble yet!!')
    
    K_,F_,Dvec,index = aModel.EliminateMatrix(lc,path)
    try:
        #sparse matrix solution         
        delta_ = sl.spsolve(K_,F_)
        delta = delta_
        f=np.zeros(len(aModel.beams)*12)
        
        #fill original displacement vector
        prev = 0
        for idx in index:
            gap=idx-prev
            if gap>0:
                delta=np.insert(delta,prev,[0]*gap)
            prev = idx + 1               
            if idx==index[-1] and idx!=len(aModel.nodes)-1:
                delta = np.insert(delta,prev, [0]*(len(aModel.nodes)*6-prev))
        delta += Dvec

        #calculate element displacement and forces
        for beam in aModel.beams.values():
            Kij_=beam.StaticCondensation()
            uij=np.zeros(12)
            fij=np.zeros(12)
            
            iend=beam.nodeI.idx
            jend=beam.nodeJ.idx

            uij[:6]=delta[iend*6:iend*6+6]
            uij[6:]=delta[jend*6:jend*6+6]
            uij = np.dot(beam.TransformMatrix(),uij)
            
            fij = np.dot(Kij_,uij) 
            
            if lc in beam.load.keys():
                fij+=beam.NodalForce(lc)

            for i in range(6):
                if beam.releaseI[i] == True:
                    fij[i] = 0
                if beam.releaseJ[i] == True:
                    fij[i + 6] = 0
            
            #beam.ID
            bid=beam.idx
            f[bid*12:bid*12+12] = fij
        
        for node in aModel.nodes.values():
            n=node.idx
            print("Disp of node "+str(n)+':')
            for i in range(n * 6,n * 6 + 6):
               print("delta[%d"%(i - n * 6) +"]=%f"%delta[i])

        for beam in aModel.beams.values():
            n=beam.idx
            print("Force of beam " +str(n)+ ':')
            for i in range(n*12,n*12+12):
                print("f[%d"%(i - n * 12)+"]=%f"%f[i])
    except Exception as e:
        print(e)
        return False
    return True
    