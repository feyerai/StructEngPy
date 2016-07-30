# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 21:57:50 2016

@author: HZJ
"""
import numpy as np
import scipy.sparse as sp

class CoordinateSystem(object):
    def __init__(self,origin, pt1, pt2):
        """
        origin: 3x1 vector
        pt1: 3x1 vector
        pt2: 3x1 vector
        """
        self.origin=origin    
        vec1 = np.array([pt1[0] - origin[0] , pt1[1] - origin[1] , pt1[2] - origin[2]])
        vec2 = np.array([pt2[0] - origin[0] , pt2[1] - origin[1] , pt2[2] - origin[2]])
        cos = np.dot(vec1, vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2)
        if  cos == 1 or cos == -1:
            raise Exception("Three points should not in a line!!")        
        self.x = vec1/np.linalg.norm(vec1)
        z = np.cross(vec1, vec2)
        self.z = z/np.linalg.norm(z)
        self.y = np.cross(self.z, self.x)
    
    def SetBy3Pts(self,origin, pt1, pt2):
        """
        origin: tuple 3
        pt1: tuple 3
        pt2: tuple 3
        """
        self.origin=origin    
        vec1 = np.array([pt1[0] - origin[0] , pt1[1] - origin[1] , pt1[2] - origin[2]])
        vec2 = np.array([pt2[0] - origin[0] , pt2[1] - origin[1] , pt2[2] - origin[2]])
        cos = np.dot(vec1, vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2)
        if  cos == 1 or cos == -1:
            raise Exception("Three points should not in a line!!")        
        self.x = vec1/np.linalg.norm(vec1)
        z = np.cross(vec1, vec2)
        self.z = z/np.linalg.norm(z)
        self.y = np.cross(self.z, self.x)
    
    def SetOrigin(self,x, y, z):
        """
        origin: tuple 3
        pt1: tuple 3
        pt2: tuple 3
        """
        self.origin = (x,y,z)
    
    def AlignWithGlobal(self):
        self.x=np.array([1,0,0])
        self.y=np.array([0,1,0])
        self.z=np.array([0,0,1])
    
    def TransformMatrix(self):
        x=self.x
        y=self.y
        z=self.z
        V=np.array([[x[0],y[0],z[0]],
                  [x[1],y[1],z[1]],
                  [x[2],y[2],z[2]]])
        return V.transpose()
        
class Material(object):
    def __init__(self,E,mu,gamma,alpha):
        """
        E\n
        mu\n
        gamma\n
        alpha
        """
        self.E = E
        self.mu = mu
        self.gamma = gamma
        self.alpha = alpha
        self.shearModulus = E / 2 / (1 + mu)

    def G(self):
        return self.shearModulus
        
class Section(object):
    def __init__(self,mat, A, J, I33, I22):
        self.material = mat
        self.A = A
        self.J = J
        self.I33 = I33
        self.I22 = I22
        
class LoadCase(object):
    def __init__(self,name,loadType):
        """
        name: name of the load case\n
        type: static,modal,spectrum,time-history
        """
        self.name=name
        self.loadType=loadType
        
class Node(object):
    def __init__(self,idx,x,y,z):
        """
        x,y,z: coordinates of nodes
        """
        self.idx=idx
        self.x=x
        self.y=y
        self.z=z
        o=[x,y,z]
        pt1=[x+1,y,z]
        pt2=[x,y+1,z]
        self.localCsys=CoordinateSystem(o,pt1,pt2)
        self.restraints=[False]*6
        self.load={}
        self.disp={}
        
    def TransformMatrix(self):
        V=self.localCsys.TransformMatrix()
        V_=np.zeros((6,6))
        V_[:3,:3]=V_[3:,3:]=V
        return V_

    def InitializeCsys(self):
        self.localCsys.AlignWithGlobal();

    def SetLoad(self,lc,load):
        """
        lc: loadcase\n
        load: a number vector indicates a nodal load.
        """
        self.load[lc]=load

    def SetDisp(self,lc,disp):
        """
        lc: loadcase\n
        load: a number vector indicates a nodal displacement.
        """
        self.disp[lc]=disp

    def SetRestraints(self,res):
        """
        res: a boolean vector indicates a nodal displacement.
        """
        self.restraints=res
        
class Beam:
    def __init__(self,idx, i, j, sec):
        """
        idx: index of the beam
        i: start node
        j: end node
        """
        self.idx=idx
        self.nodeI=i
        self.nodeJ=j
        self.load={}
        self.releaseI=[False]*6
        self.releaseJ=[False]*6
        self.section=sec
        tol = 1E-6
        #Initialize local CSys
        o = [ self.nodeI.x, self.nodeI.y, self.nodeI.z ]
        pt1 = [ self.nodeJ.x, self.nodeJ.y, self.nodeJ.z ]
        pt2 = [ self.nodeI.x, self.nodeI.y, self.nodeI.z ]
        if abs(self.nodeI.x - self.nodeJ.x) < tol and abs(self.nodeI.y - self.nodeJ.y) < tol:
            pt2[0] += 1
        else:
            pt2[2] += 1
        self.localCsys = CoordinateSystem(o, pt1, pt2)

        #Initialize local stiffness matrix
        l = self.Length()
        E = self.section.material.E
        A = self.section.A
        J = self.section.J
        G = self.section.material.G()
        I2 = self.section.I22
        I3 = self.section.I33
        rho = self.section.material.gamma

        self.Kij = np.zeros((12, 12))
        self.Mij = np.zeros((12, 12))

        #form the stiffness matrix:
        self.Kij[0, 0]=E*A / l
        self.Kij[0, 6]=self.Kij[6, 0]=-E*A / l

        self.Kij[1, 1]=12 * E*I3 / l / l / l
        self.Kij[1, 5]=self.Kij[5, 1]=6 * E*I3 / l / l
        self.Kij[1, 7]=self.Kij[7, 1]=-12 * E*I3 / l / l / l
        self.Kij[1, 11]=self.Kij[11, 1]=6 * E*I3 / l / l

        self.Kij[2, 2]=12 * E*I2 / l / l / l
        self.Kij[2, 4]=self.Kij[4, 2]=-6 * E*I2 / l / l
        self.Kij[2, 8]=self.Kij[8, 2]=-12 * E*I2 / l / l / l
        self.Kij[2, 10]=self.Kij[10, 2]=-6 * E*I2 / l / l

        self.Kij[3, 3]=G*J / l
        self.Kij[3, 9]=self.Kij[9, 3]=-G*J / l

        self.Kij[4, 4]=4 * E*I2 / l
        self.Kij[4, 8]=self.Kij[8, 4]=6 * E*I2 / l / l
        self.Kij[4, 10]=self.Kij[10, 4]=2 * E*I2 / l

        self.Kij[5, 5]=4 * E*I3 / l
        self.Kij[5, 7]=self.Kij[7, 5]=-6 * E*I3 / l / l
        self.Kij[5, 11]=self.Kij[11, 5]=2 * E*I3 / l

        self.Kij[6, 6]=E*A / l

        self.Kij[7, 7]=12 * E*I3 / l / l / l
        self.Kij[7, 11]=self.Kij[11, 7]=-6 * E*I3 / l / l

        self.Kij[8, 8]=12 * E*I2 / l / l / l
        self.Kij[8, 10]=self.Kij[10, 8]=6 * E*I2 / l / l

        self.Kij[9, 9]=G*J / l

        self.Kij[10, 10]=4 * E*I2 / l

        self.Kij[11, 11]=4 * E*I3 / l

        #form mass matrix    
        ##Coordinated mass matrix
        #Mij[0, 0]=140
        #Mij[0, 6]=70

        #Mij[1, 1]=156
        #Mij[1, 5]=Mij[5, 1]=22 * l
        #Mij[1, 7]=Mij[7, 1]=54
        #Mij[1, 11]=Mij[11, 1]=-13 * l

        #Mij[2, 2]=156
        #Mij[2, 4]=Mij[4, 2]=-22 * l
        #Mij[2, 8]=Mij[8, 2]=54
        #Mij[2, 10]=Mij[10, 2]=13 * l

        #Mij[3, 3]=140 * J / A
        #Mij[3, 9]=Mij[9, 3]=70 * J / A

        #Mij[4, 4]=4 * l *l
        #Mij[4, 8]=Mij[8, 4]=-13 * l
        #Mij[4, 10]=Mij[10, 4]=-3 * l*l

        #Mij[5, 5]=4 * l*l
        #Mij[5, 7]=Mij[7, 5]=13 * l
        #Mij[5, 11]=Mij[11, 5]=-3 * l*l

        #Mij[6, 6]=140

        #Mij[7, 7]=156
        #Mij[7, 11]=Mij[11, 7]=-22 * l

        #Mij[8, 8]=156
        #Mij[8, 10]=Mij[10, 8]=22 * l

        #Mij[9, 9]=140 * J / A

        #Mij[10, 10]=4 * l*l

        #Mij[11, 11]=4 * l*l

        #Mij*= (rho*A*l / 420)

        #Concentrated mass matrix
        for i in range(12):
            self.Mij[i, i]=1
        self.Mij*=rho*A*l/2

    def InitializeCsys(self):
        nodeI=self.nodeI
        nodeJ=self.nodeJ
        o = np.array([nodeI.x, nodeI.y, nodeI.z])
        pt1 = np.array([nodeJ.x, nodeJ.y, nodeJ.z])
        pt2 = np.array([0,0,0])
        if self.nodeI.x != self.nodeJ.x and self.nodeI.y != self.nodeJ.y:
            pt2[2] = 1
        else:
            pt2[0] = 1
        self.localCsys.SetBy3Pts(o, pt1, pt2)

    def Length(self):
        nodeI=self.nodeI
        nodeJ=self.nodeJ
        return np.sqrt((nodeI.x - nodeJ.x)*(nodeI.x - nodeJ.x) + (nodeI.y - nodeJ.y)*(nodeI.y - nodeJ.y) + (nodeI.z - nodeJ.z)*(nodeI.z - nodeJ.z))

    #vec NodalForceI()
    #{
    #    l = self.Length()
    #    #recheck!!!!!!!!!!!!
    #    vec v(6, fill::zeros)
    #    v(0) = (loadI[0] + loadJ[0]) * l / 2#P
    #    v(1) = (loadI[1] * 7 / 20 + loadJ[1] * 3 / 20) * l#V2
    #    v(2) = (loadI[2] * 7 / 20 + loadJ[2] * 3 / 20) * l#V3
    #    v(3) = loadI[3] - loadJ[3]#T
    #    v(4) = (loadI[2] / 20 + loadJ[2] / 30) * l * l + loadI[4]#M22
    #    v(5) = (loadI[1] / 20 + loadJ[1] / 30) * l * l + loadI[5]#M33
    #
    #    return v
    #    #l = self.Length()
    #    #mat bij(6, 6,fill::zeros)
    #    #for (int i = 0 i < 6 i++)
    #    #    bij(i, i) = -1
    #    #bij(1, 5) = bij(2, 4) = 1/l
    #    #bij(5, 1) = bij(4, 2) = l
    #    ###############TEST##############
    #    #mat k = NodalForceJ()
    #    ##cout.setf(ios::scientific)
    #    #for (int i = 0 i < k.n_rows i++)
    #    #{
    #    #    for (int j = 0 j < k.n_cols j++)
    #    #    {
    #    #        cout.width(8)
    #    #        cout.precision(4)
    #    #        cout.setf(ios::right)
    #    #        cout << k[i, j] << "  "
    #    #    }
    #    #    cout << endl
    #    #}
    #    #cout << endl
    #    ###############TEST##############
    #    ###############TEST##############
    #    #k = bij
    #    ##cout.setf(ios::scientific)
    #    #for (int i = 0 i < k.n_rows i++)
    #    #{
    #    #    for (int j = 0 j < k.n_cols j++)
    #    #    {
    #    #        cout.width(8)
    #    #        cout.precision(4)
    #    #        cout.setf(ios::right)
    #    #        cout << k[i, j] << "  "
    #    #    }
    #    #    cout << endl
    #    #}
    #    #cout << endl
    #    ###############TEST##############
    #    ###############TEST##############
    #    # k = bij*NodalForceJ()
    #    ##cout.setf(ios::scientific)
    #    #for (int i = 0 i < k.n_rows i++)
    #    #{
    #    #    for (int j = 0 j < k.n_cols j++)
    #    #    {
    #    #        cout.width(8)
    #    #        cout.precision(4)
    #    #        cout.setf(ios::right)
    #    #        cout << k[i, j] << "  "
    #    #    }
    #    #    cout << endl
    #    #}
    #    #cout << endl
    #    ###############TEST##############
    #    #return bij*NodalForceJ()
    #}

    def NodalForce(self,lc):
        """
        returns a 12x1 element nodal force vector
        """
        l = self.Length()
        loadI=self.load[lc][:6]
        loadJ=self.load[lc][6:]
        #recheck!!!!!!!!!!!!
        #i
        v=np.zeros(12)
        v[0]=(loadI[0] + loadJ[0]) * l / 2#P
        v[1]=(loadI[1] * 7 / 20 + loadJ[1] * 3 / 20) * l#V2
        v[2]=(loadI[2] * 7 / 20 + loadJ[2] * 3 / 20) * l#V3
        v[3]=loadI[3] - loadJ[3]#T
        v[4]=(loadI[2] / 20 + loadJ[2] / 30) * l * l + loadI[4]#M22
        v[5]=(loadI[1] / 20 + loadJ[1] / 30) * l * l + loadI[5]#M33
        #j
        v[6]=(loadJ[0] + loadI[0]) * l / 2#P
        v[7]=(loadJ[1] * 7 / 20 + loadI[1] * 3 / 20) * l#V2
        v[8]=(loadJ[2] * 7 / 20 + loadI[2] * 3 / 20) * l#V3
        v[9] = loadJ[3] - loadI[3]#T
        v[10] = -(loadJ[2] / 20 + loadI[2] / 30) * l * l + loadJ[4]#M22
        v[11] = -(loadJ[1] / 20 + loadI[1] / 30) * l * l + loadJ[5]#M33
        return v

    def LocalStiffnessMatrix(self):
        return self.Kij

    def LocalMassMatrix(self):
        return self.Mij

    def StaticCondensation(self, mass=False):
        """
        return:
        kij_: 12x12 matrix
        mij_: 12x12 matrix, if mass==True
        """
        if mass==False:
            kij=self.Kij
            kij_ = kij
    
            for n in range(6):
                if self.releaseI[n] == True:
                    for i in range(12):
                        for j in range(12):
                            kij_[i, j] = kij[i, j] - kij[i, n]* kij[n, j] / kij[n, n]
                if self.releaseJ[n] == True:
                    for i in range(12):
                        for j in range(12):
                            kij_[i, j] = kij[i, j] - kij[i, n + 6]* kij[n + 6, j] / kij[n + 6, n + 6]
            return kij_
        else:
            kij=self.Kij
            mij=self.Mij
            kij_ = kij.copy()
            mij_ = mij.copy()
    
            for n in range(0,6):
                if self.releaseI[n] == True:
                    for i in range(12):
                        for j in range(12):
                            kij_[i, j] = kij[i, j] - kij[i, n]* kij[n, j] / kij[n, n]
                            mij_[i, j] = mij[i, j] - mij[i, n]* mij[n, j] / mij[n, n]
                if self.releaseJ[n] == True:
                    for i in range(12):
                        for j in range(12):
                            kij_[i, j] = kij[i, j] - kij[i, n + 6]* kij[n + 6, j] / kij[n + 6, n + 6]
                            mij_[i, j] = mij[i, j] - mij[i, n + 6]* mij[n + 6, j] / mij[n + 6, n + 6]
            return kij_, mij_
            
    def LoadCondensation(self,lc):
        """
        lc: loadcase
        return???refer????
        """
        kij=self.Kij
        rij=self.NodalForce(lc)
        kij_ = kij.copy()
        rij_ = rij.copy()

        for n in range(6):
            if self.releaseI[n] == True:
                for i in range(12):
                    for j in range(12):
                        kij_[i, j] = kij[i, j] - kij[i, n]* kij[n, j] / kij[n, n]
                    rij_[i] = rij[i] - rij[n] * kij[n, i] / kij[n, n]
            if self.releaseJ[n] == True:
                for i in range(12):
                    for j in range(12):
                        kij_[i, j] = kij[i, j] - kij[i, n + 6]* kij[n + 6, j] / kij[n + 6, n + 6]
                    rij_[i] = rij[i] - rij[n + 6] * kij[n + 6, i] / kij[n + 6, n + 6]
        return rij_

    def CalculateElmForce(self,uij,fij):
        """
        uij,fij: 12x1 sparse vector
        """
        fij = np.zeros(12)
        Kij = sp.csc_matrix(12, 12)
        rij = sp.csc_matrix(12,1)
        self.StaticCondensation(Kij, rij)
        fij = Kij * uij + self.NodalForce()
        return fij

    def TransformMatrix(self):
        T=np.zeros((12,12))
        V=self.localCsys.TransformMatrix()
        T[:3,:3] =T[3:6,3:6]=T[6:9,6:9]=T[9:,9:]= V
        return T

    def SetLoadDistributed(self,lc,q):
        """
        lc: loadcase\n
        q: 12x1 vector represent distributed forces
        """
        self.load[lc]=q
        