# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 21:32:57 2016

@author: HZJ
"""
import sqlite3
import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as sl
from scipy import linalg
import pandas as pd
import Material,Section,Node,Beam
import DataManager

class Model:
    def __init__(self,db):
        self.data=DataManager.ModelData.ModelData()
        self.database=db
        
    def Save(self):
        self.data.Save(self.database)
        
    def Run(self,loadcases):
        """
        loadcases: load cases set to run.
        """
        model=AnalysisModel(self.data)
        model.Assemble()
        for lc in loadcases:
            model.AssembleLoad(lc)
            model.SolveLinear(lc)
        
    def Test(self):
        #*************************************************Modeling**************************************************/
        #Q345
        materials=[]
        sections=[]
        nodes=[]
        beams=[]
        
        materials.append(Material.Material(2.000E11, 0.3, 7849.0474, 1.17e-5))
        #H200x150x8x10
        sections.append(Section.Section(materials[0], 4.800E-3, 1.537E-7, 3.196E-5, 5.640E-6))
        nodes.append(Node.Node(0, 0, 0, 0))
        nodes.append(Node.Node(1, 5, 0, 0))
        nodes.append(Node.Node(2, 10,0, 0))
        for i in range(len(nodes)-1):
            beams.append(Beam.Beam(i, nodes[i], nodes[i+1], sections[0]))
        qi=(0,0,-10,0,0,0)
        qj=(0,0,-10,0,0,0)
        beams[0].SetLoadDistributed(qi, qj)
        res1 = [True]*6
        res2 = [False]*6
        res3 = [False,True,True,True,True,True]
        nodes[0].SetRestraints(res1)
        nodes[1].SetRestraints(res2)
        #*************************************************Modeling**************************************************/
        ID=0
        ns=[]
        for node in nodes:
            ns.append([ID,node.x,node.y,node.z])
            ID+=1
            
        md=self.data
        #Define material
        DataManager.Definition.Properties.Materials.AddQuick(md,'GB','Q345')
        #Define section
        DataManager.Definition.Properties.BeamSections.AddQuick(md,'Q345','H400x200x12x14')
        #Define load case
        DataManager.Definition.LoadCases.StaticLinear.SetCase(md,'SD')
        #Add nodes
        DataManager.Modeling.Nodes.AddCartesian(md,ns)
        #Add beams
        DataManager.Modeling.Beams.AddBeams(md,[(0,0,1,'H400x200x12x14')])
        DataManager.Modeling.Beams.AddBeams(md,[(1,1,2,'H400x200x12x14')])
        #Set loads
        DataManager.Modeling.Nodes.SetLoadForce(md,1,'SD',[0,0,-10,0,0,0])
        #Set restraints
        DataManager.Modeling.Nodes.SetRestraints(md,0,res2)
        DataManager.Modeling.Nodes.SetRestraints(md,2,res2)
#        self.Save()
        aModel=AnalysisModel(self.data)
        path='F:\\Test\\'
        aModel.Assemble(path)
        aModel.AssembleLoad(path)

        
class AnalysisModel(object):
    """
    Assembler handles the data read from DataManager and produces nodes and elements.
    """        
    def __init__(self,md): 
        self.materials=self.GetMaterials(md)                
        self.beamSections=self.GetBeamSections(md)
        self.nodes=self.GetNodes(md)
        self.beams=self.GetBeams(md)
        self.quads=self.GetQuads(md)
        
#        for idx,row in self.nodedf.iterrows():
#            node=Node.Node(row['X'],row['Y'],row['Z'])
#            hid=row['HID']
#            self.nodes[hid]=node
#            
#        for idx,row in beamdf.iterrows():
#            nodeI=self.node[row['NodeI']]
#            nodeJ=self.node[row['NodeJ']]
#            section=row['Section']
#            beam=Beam.Beam(nodeI,nodeJ,)
#            hid=row['HID']
#            self.beams[hid]=beam
    def GetMaterials(self,md):
        """
        returns a dictionary of node objects
        """
        mg=md.dataFrames['MaterialPropertiesGeneral']
        ms=md.dataFrames['MaterialPropertiesSteel']
        mb=md.dataFrames['MaterialPropertiesBasicMechanical']
        df=pd.DataFrame.join(mg,mb)
        df=pd.DataFrame.join(df,ms)
        materials={}
        for idx,row in df.iterrows():
            material=Material.Material(row['E1'],row['G12'],row['UnitWeight'],row['A1'])
            materials[idx]=material
        return materials
        
    def GetBeamSections(self,md):
        """
        returns a dictionary of node objects, should be executed after GetMaterials()
        """
        df=md.dataFrames['BeamPropGeneral']
        beamSections={}
        for idx,row in df.iterrows():
            beamSection=Section.Section(self.materials[row['Material']],row['Area'],row['TorsConst'],row['I33'],row['I22'])
            beamSections[idx]=beamSection
        return beamSections
        
    def GetNodes(self,md):
        """
        returns a dictionary of node objects
        """
        coor=md.dataFrames['NodeCoordinates']
        rest=md.dataFrames['NodeRestraintAssignment']
        nodes={}
        HID=0
        for idx,row in coor.iterrows():
            node=Node.Node(HID,row['X'],row['Y'],row['Z'])
            nodes[idx]=node
            HID+=1           
        for idx,row in rest.iterrows():
            res=[]
            for i in row:
                if i:
                    res.append(False)
                else:
                    res.append(True)
            nodes[idx].SetRestraints(res)
        return nodes
        
    def GetBeams(self,md):
        """
        returns a dictionary of node objects, should be executed after GetNodes(), GetMaterials() and GetBeamSections()
        """
        beams={}
        bc=md.dataFrames['ConnectivityBeam']
        bs=md.dataFrames['BeamSectionAssignments']
        df=pd.DataFrame.join(bc,bs)
        HID=0
        for idx,row in df.iterrows():
            beam=Beam.Beam(HID,self.nodes[row['NodeI']],self.nodes[row['NodeJ']],self.beamSections[row['Section']])
            beams[idx]=beam
            HID+=1
        return beams
        
    def GetQuads(self,md):
        return False

    def Update(self):
        """
        Update model for multi-step or non-linear analysis
        """
        return False
            
    def Assemble(self,path):
        """
        Assemble matrix
        Writing the matrix to the disk
        """        
        #*************************************************Modeling**************************************************/
        #use dictionary to improve random access performance
        n=len(self.nodes.keys())
        # Dynamic space allocate
        Kmat = np.zeros((n*6, n*6))
        Mmat = np.zeros((n*6, n*6))
        Dvec = np.zeros(n*6)
            
        #Beam load and displacement, and reset the index
        for key,beam in self.beams.items():
            i = beam.nodeI.idx
            j = beam.nodeJ.idx
            T=np.matrix(beam.TransformMatrix())
            Tt = T.transpose()

            #Transform matrix
            Vl=np.matrix(beam.localCsys.TransformMatrix())
            V=np.zeros((6, 6))
            V[:3,:3] =V[3:,3:]= Vl
            Vt = V.T

            #Static condensation to consider releases
            Kij=np.zeros((12, 12))
            Mij=np.zeros((12, 12))
            rij=np.zeros(12)
            Kij, rij, Mij = beam.StaticCondensation(Kij, rij, Mij)

            #Assemble Total Stiffness Matrix
            Ke = np.dot(np.dot(Tt,Kij),T)
            Keii = Ke[:6,:6]
            Keij = Ke[:6,6:]
            Keji = Ke[6:,:6]
            Kejj = Ke[6:,6:]
            Kmat[i*6:i*6+6, i*6:i*6+6] += Keii
            Kmat[i*6:i*6+6, j*6:j*6+6] += Keij
            Kmat[j*6:j*6+6, i*6:i*6+6] += Keji
            Kmat[j*6:j*6+6, j*6:j*6+6] += Kejj

            #Assembel Mass Matrix        
            Me = np.dot(np.dot(Tt,Mij),T)
            Meii = Me[:6,:6]
            Meij = Me[:6,6:]
            Meji = Me[6:,:6]
            Mejj = Me[6:,6:]
            Mmat[i*6:i*6+6, i*6:i*6+6] += Meii
            Mmat[i*6:i*6+6, j*6:j*6+6] += Meij
            Mmat[j*6:j*6+6, i*6:i*6+6] += Meji
            Mmat[j*6:j*6+6, j*6:j*6+6] += Mejj
        
        for key,node in self.nodes.items():
            nid=node.idx
            for i in range(6):
                if node.disp[i] != 0:
                    Dvec[nid * 6 + i] = node.disp[i]
                    
        scipy.io.mmwrite(path+'K.mtx',sp.coo_matrix(Kmat))
        scipy.io.mmwrite(path+'M.mtx',sp.coo_matrix(Mmat))
        scipy.io.mmwrite(path+'D.mtx',sp.coo_matrix(Dvec))
        
    def AssembleLoad(self,path,lc=''):
        """
        Assemble load vector for specified load case
        lc: loadcase
        """
        n=len(self.nodes.keys())
        Fvec = np.zeros(n*6)
        
        for key,node in self.nodes.items():
            nid=node.idx
            load = np.array(node.load)
            Fvec[nid * 6: nid * 6 + 6] = np.dot(node.TransformMatrix().T,load)
                 
        for key,beam in self.beams.items():
            #Transform matrix
            Vl=np.matrix(beam.localCsys.TransformMatrix())
            V=np.zeros((6, 6))
            V[:3,:3] =V[3:,3:]= Vl
            Vt = V.T
            #Static condensation to consider releases
            Kij=np.zeros((12, 12))
            Mij=np.zeros((12, 12))
            rij=np.zeros(12)
            Kij, rij, Mij = beam.StaticCondensation(Kij, rij, Mij)
            #Assemble nodal force vector
            i = beam.nodeI.idx
            j = beam.nodeJ.idx
            Fvec[i*6:i*6+6] += np.dot(Vt,rij[:6])
            Fvec[j*6:j*6+6] += np.dot(Vt,rij[6:])           
        scipy.io.mmwrite(path+'F.mtx',sp.coo_matrix(Fvec))

    def isAssembled(self):
#        if self.Kmat == None:
#            return False
        return True

    def K(self):
        return self.Kmat

    def F(self):
        return self.Fvec

    def EliminateMatrix(self,path,mass=False):
        """
        return 
        K_bar: sparse matrix
        F_bar: sparse matrix
        M_bar: sparse matrix
        index: vector
        """
        Kmat=scipy.io.mmread(path+'K.mtx')
        Mmat=scipy.io.mmread(path+'M.mtx')
        Fvec=scipy.io.mmread(path+'F.mtx').toarray()[0]
        Dvec=scipy.io.mmread(path+'D.mtx').toarray()[0]
        if mass==False:
            k = Kmat.todense()
            f = Fvec
            Id=np.arange(len(f))
            nRemoved=0
            i=0
            for node in self.nodes:
                for j in range(6):
                    if node.restraints[j] == True or node.disp[j] != 0:
                        k=np.delete(k,i*6+j-nRemoved,axis=0)
                        k=np.delete(k,i*6+j-nRemoved,axis=1)
                        f=np.delete(f,i*6+j-nRemoved)
                        Id=np.delete(Id,i*6+j-nRemoved)
                        nRemoved+=1
                i+=1
            K_ = k
            F_ = f
            index = Id
            return K_,F_,Dvec,index
        else:
            k = Kmat
            m = Mmat
            f = Fvec
            Id=np.arange(len(f))
            nRemoved = 0
            i=0
            for node in self.nodes:
                for j in range(6):
                    if node.restraints[j] == True or node.disp[j]!=0:
                        k=np.delete(k,i*6+j-nRemoved,axis=0)
                        k=np.delete(k,i*6+j-nRemoved,axis=1)
                        f=np.delete(f,i*6+j-nRemoved)
                        m=np.delete(m,i*6+j-nRemoved,axis=0)
                        m=np.delete(m,i*6+j-nRemoved,axis=1)
                        Id=np.delete(Id,i*6+j-nRemoved)
                        nRemoved+=1
                i+=1
            K_ = k
            M_ = m
            F_ = f
            index = Id
            return K_,M_,F_,Dvec,index

    def SolveLinear(self,path):
        if not self.isAssembled():
            raise Exception('Not assemble yet!!')
        
        K_bar,F_bar,Dvec,index = self.EliminateMatrix(path)
        try:
            #sparse matrix solution         
            delta_bar = sl.spsolve(K_bar,F_bar)
            print('HERE!')
            delta = delta_bar
            f=np.zeros(len(self.beams)*12)
           
            #fill original displacement vector
            prev = 0
            for idx in index:
                gap=idx-prev
                if gap>0:
                    delta=np.insert(delta,prev,[0]*gap)
                prev = idx + 1               
                if idx==index[-1] and idx!=len(self.nodes)-1:
                    delta = np.insert(delta,prev, [0]*(len(self.nodes)*6-prev))
            delta += Dvec

            #calculate element displacement and forces
            for beam in self.beams:
                Kij_bar=np.zeros((12, 12))
                rij_bar=np.zeros((12,1))
                Kij_bar,rij_bar=beam.StaticCondensation(Kij_bar, rij_bar)
                uij=np.zeros(12)
                fij=np.zeros(12)
                
                i=0
                for node in self.nodes:
                    if node is beam.nodeI:
                        iend=i
                    i+=1
                i=0
                for node in self.nodes:
                    if node is beam.nodeJ:
                        jend=i
                    i+=1
                
                uij[:6]=delta[iend*6:iend*6+6]
                uij[6:]=delta[jend*6:jend*6+6]
                uij = np.dot(beam.TransformMatrix(),uij)
                
                fij = np.dot(Kij_bar,uij) + beam.NodalForce()
                for i in range(6):
                    if beam.releaseI[i] == True:
                        fij[i] = 0
                    if beam.releaseJ[i] == True:
                        fij[i + 6] = 0
                #beam.ID
                i=0
                for b in self.beams:
                    if beam is b:
                        bid=i
                    i+=1
                f[bid*12:bid*12+12] = fij
            for n in range(len(self.nodes)):
                print("Disp of node "+str(n)+':')
                for i in range(n * 6,n * 6 + 6):
                   print("delta[%d"%(i - n * 6) +"]=%f"%delta[i])

            for n in range(len(self.beams)):
                print("Force of beam " +str(n)+ ':')
                for i in range(n*12,n*12+12):
                    print("f[%d"%(i - n * 12)+"]=%f"%f[i])
        except Exception as e:
            print(e)
            return False
        return True

    def SolveModal(self,path,k):
        if not self.isAssembled():
            raise Exception('Not assemble yet!!')           
        K_bar,M_bar,F_bar,index,Dvec = self.EliminateMatrix(path,True)

        try:
            #general eigen solution, should be optimized later!!
            A=sp.csr_matrix(K_bar)
            B=sp.csr_matrix(M_bar)
            
            omega2s,mode = sl.eigsh(A,k,B)
        
            for omega2 in omega2s:
                print(2*3.14/np.sqrt(omega2))

            #extract vibration mode
            delta_bar = np.linalg.norm(mode)/np.linalg.norm(mode)
            delta=delta_bar
            f=np.zeros(len(self.beams)*12)
            
            print('HERE!!!!')
            #fill original displacement vector
            prev = 0
            for idx in index:
                gap=idx-prev
                if gap>0:
                    delta=np.insert(delta,prev,[0]*gap)
                prev = idx + 1               
                if idx==index[-1] and idx!=len(self.nodes)-1:
                    delta = np.insert(delta,prev, [0]*(len(self.nodes)*6-prev))
            delta += Dvec

            #calculate element displacement and forces
            for beam in self.beams:
                Kij_bar=np.zeros((12, 12))
                rij_bar=np.zeros(12)
                Kij_bar,rij_bar=beam.StaticCondensation(Kij_bar,rij_bar)
                uij=np.zeros(12)
                fij=np.zeros(12)
                
                i=0
                for node in self.nodes:
                    if node is beam.nodeI:
                        iend=i
                    i+=1
                i=0
                for node in self.nodes:
                    if node is beam.nodeJ:
                        jend=i
                    i+=1
                
                uij[:6]=delta[iend*6:iend*6+6]
                uij[6:]=delta[jend*6:jend*6+6]
                uij = np.dot(beam.TransformMatrix(),uij)

                fij = np.dot(Kij_bar,uij) + beam.NodalForce()
                for i in range(6):
                    if beam.releaseI[i] == True:
                        fij[i] = 0
                    if beam.releaseJ[i] == True:
                        fij[i + 6] = 0
                
                #beam.ID
                i=0
                for b in self.beams:
                    if beam is b:
                        bid=i
                    i+=1  
            
                f[bid * 12,bid*12+12]=fij
            for n in range(len(self.nodes)):
                print("Disp of node " +str(n)+ ':')
                for i in range(n*6, n*6+6):
                    print("delta["+str(i-n*6)+"]="+str(delta[i]))
            for n in range(len(self.beams)):
                print( "Force of beam "+str(n)+':')
                for i in range(n*12, n*12 + 12):
                    print("f[" +str(i-n*12)+"]=" +str(f[i]))
        except Exception as e:
            print(str(e))
            return False
        return True

    def SetMass(self):
        for beam in self.beams:
            beam.nodeI.mass += beam.section.A*beam.section.material.gamma*beam.Length() / 2
            beam.nodeJ.mass += beam.section.A*beam.section.material.gamma*beam.Length() / 2
        return False        
        
        
        
if __name__=='__main__':     
    file='F:\\Test\\Test.sqlite'
    path='F:\\Test'
    m=Model(file)
    m.Test()
#    m.Assemble(path)
#    m.SolveLinear(path)
#    m.SolveModal(path,k=1)
