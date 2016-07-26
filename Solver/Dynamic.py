# -*- coding: utf-8 -*-
"""
Created on Sun May 29 14:01:08 2016

@author: HZJ
"""

import numpy as np
import pandas as pd
import sys
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.fftpack import fft, ifft

class SDOFSystem():
    def __init__(self,m,k,zeta):
        """
        m: mass of the system\n
        k: stiffness of the system\n
        zeta: damping ratio
        """
        self.m=m
        self.k=k
        self.zeta=zeta
        
    def omega_n(self):
        """
        frequency without damping
        """
        return np.sqrt(self.k/self.m)
        
    def omega_d(self):
        """
        frequency with low damping
        """
        return self.omega_n()*np.sqrt(1-self.zeta**2)  
        
    def Duhamel(self,t,p,u0=0,v0=0):
        """
        Convolutional Integration to solve time history excitation.
        """
        m=self.m
        zeta=self.zeta
        omega_d=self.omega_d()
        omega_n=self.omega_n()
        #discreted convolution
        convol=0
        for tau in range(len(t)-1):
            dtau=t[tau+1]-t[tau]
            convol+=p[tau]*np.exp(-zeta*omega_n*(t[-1]-t[tau]))*np.sin(omega_d*(t[-1]-t[tau]))*dtau
        convol*=(1/m/omega_d)
        #effect of initial state
        ini=u0*np.exp(-zeta*omega_n*t[-1])*np.cos(omega_d*t[-1])+1/omega_d*(v0+zeta*omega_n*u0)*np.exp(-zeta*omega_n*t[-1])*np.sin(omega_d*t[-1])       
        return convol+ini
        
    def InterpExcite(self,t,p,u0=0,v0=0):
        """
        Solve time history excitation with inpterpolate excite method.\n
        t: time\n
        p: time-relative exciting force\n
        optional:\n
        u0: initial displacement.\n
        v0: initial velocity.\n
        """
        k=self.k
        zeta=self.zeta
        omega_d=self.omega_d()
        omega_n=self.omega_n()
        beta=zeta*omega_n
        h=0
        u=[u0]
        v=[v0]
        for i in range(len(t)-1):
            if t[i+1]-t[i]!=h:
                h=t[i+1]-t[i]
                A=1/(k*omega_d*h)*(np.exp(-beta*h)*(((omega_d**2-beta**2)/omega_n**2-beta*h)*np.sin(omega_d*h)-(2*omega_d*beta/omega_n**2+omega_d*h)*np.cos(omega_d*h))+2*beta*omega_d/omega_n**2)
                B=1/(k*omega_d*h)*(np.exp(-beta*h)*(-(omega_d**2-beta**2)/omega_n**2*np.sin(omega_d*h)+(2*omega_d*beta/omega_n**2)*np.cos(omega_d*h))+omega_d*h-2*beta*omega_d/omega_n**2)
                C=np.exp(-beta*h)*(np.cos(omega_d*h)+beta/omega_d*np.sin(omega_d*h))
                D=1/omega_d*np.exp(-beta*h)*np.sin(omega_d*h)        
                A1=1/(k*omega_d*h)*(np.exp(-beta*h)*((beta+omega_n**2*h)*np.sin(omega_d*h)+omega_d*np.cos(omega_d*h))-omega_d)        
                B1=1/(k*omega_d*h)*(-np.exp(-beta*h)*(beta*np.sin(omega_d*h)+omega_d*np.cos(omega_d*h))+omega_d)        
                C1=-(omega_n**2/omega_d)*np.exp(-beta*h)*np.sin(omega_d*h)
                D1=np.exp(-beta*h)*(np.cos(omega_d*h)-beta/omega_d*np.sin(omega_d*h))
            u.append(A*p[i]+B*p[i]+C*u[i]+D*v[i])
            v.append(A1*p[i]+B1*p[i+1]+C1*u[i]+D1*v[i])
        print(u)    
        df=pd.DataFrame({'t':t,'u':u,'v':v})
        return df
    
    def AverageAcceleration(self,t,p,u0=0,v0=0):
        """
        Solve time history excitation with average accerleration method.\n
        t: time\n
        p: time-relative exciting force\n
        optional:\n
        u0: initial displacement.\n
        v0: initial velocity.\n
        """
        m=self.m
        k=self.k
        zeta=self.zeta
        c=zeta*2*np.sqrt(k*m)
        dt=0
        u=[u0]
        v=[v0]
        a=[1/m*(p[0]-c*v[0]-k*u[0])]
        for i in range(len(t)-1):
            if t[i+1]-t[i]!=dt:
                dt=t[i+1]-t[i]
                k1=k+2*c/dt+4*m/dt**2
            dp=p[i+1]-p[i]
            dp1=dp+(4*m/dt+2*c)*v[i]+2*m*a[i]
            du=dp1/k1
            da=(4/dt**2)*(du-v[i]*dt)-2*a[i]
            dv=(2/dt)*du-2*v[i]
            u.append(u[i]+du)
            v.append(v[i]+dv)
            a.append(a[i]+da)
        df=pd.DataFrame({'t':t,'u':u,'v':v})
        return df

class MDOFSystem():
    def __init__(self,M,K,C):
        """
        M: Mass matrix\n
        K: Stiffness matrix\n
        C: Damping matrix\n
        """
        self.M=M
        self.K=K
        self.C=C
        
    def EigenMode(self,n):
        """
        Solve the eigen mode of the MDOF system\n
        n: number of modes to extract.
        """
#        w2,mode=eigsh(A=self.K,M=self.M,k=n,which='LM')
        w2,mode=eigh(self.K,self.M)
        freq=[]
        w=[]
        T=[]
        for i in range(n):
            w.append(np.sqrt(w2[i]))
            T.append(2*np.pi/w[-1])            
#            if min(mode[i]/max(mode[i]))<-1:
#                mode[i]/=min(mode[i])
#            else:
#                mode[i]/=max(mode[i])
            freq.append(1/T[-1])
        return w,freq,T,mode
        
    def RizMode(self,n,F):
        """
        Solve the Riz mode of the MDOF system\n
        n: number of modes to extract\n
        F: spacial load pattern
        """
        
        return False   
    
    def SpectrumAnalysis(self,n,spec):
        """
        sepctrum analysis\n
        n: number of modes to use\n
        spec: a list of tuples (period,acceleration response)
        """
        freq,mode=self.EigenMode(n)
        M_=np.dot(mode.T,self.M)
        M_=np.dot(M_,mode)
        K_=np.dot(mode.T,self.K)
        K_=np.dot(K_,mode)
        C_=np.dot(mode.T,self.C)
        C_=np.dot(C_,mode)
        d_=[]
        for (m_,k_,c_) in zip(M_.diag(),K_.diag(),C_.diag()):
            sdof=SDOFSystem(m_,k_)
            T=sdof.omega_d()
            d_.append(np.interp(T,spec[0],spec[1]*m_))
        d=np.dot(d_,mode)
        #CQC
            
        
        
        return d
        
    def ModalDecomposition(self,T,F,u0,v0,a0,xi):
        """
        Solve time-history problems with modal decomposition method.\n
        u0,v0,a0: initial state.\n
        T: time list with uniform interval.\n
        F: list of time-dependent force vectors.\n
        xi: modal damping ratio
        """
        wVec,fVec,TVec,modeMat=self.EigenMode(self.K.shape[0])
        damp=np.diag(2*xi*wVec)
        y_mat=np.dot(modeMat.T,np.array([u0,v0,a0]))
        dt=T[1]-T[0]        
        #iterative solve        
        u=[]
        v=[]
        a=[]
        for w,mode,y0 in zip(wVec,modeMat,y_mat):    
            #construct load vector R
            R=[]
            R0=np.dot(mode.T,F.T)
            R1=[]
            R2=[]
            R3=[]
            R1=(R0[1]-R0[0])/dt
            for i in range(len(T)-1):
                R1=np.append(R1,(R0[i+1]-R0[i])/dt)            
            R1=np.append(R1,[-1])
            for i in range(len(T)-1):
                R2.append(6/dt/dt*(R0[i]-R0[i+1])+2/dt*(R1[i+1]+2*R1[i]))
            R2.append(R2[-1])
            for i in range(len(T)-1):
                R3.append((R2[i+1]-R2[i])/dt) 
            R3.append(R3[-1])
            for i in range(len(T)):
                t=T[i]
                R.append(R0[i]+t*R1[i]+t**2/2*R2[i]+t**3/6*R3[i])
#            p=mode.T*f #mode paticipation
            wd=w*np.sqrt(1-xi**2)        
            w_=w*xi
            xi_=xi/np.sqrt(1-xi**2)
            a0=2*xi*w
            a1=wd**2-w_**2
            a2=2*w_*wd
            S0=np.exp(-xi*w*dt)*np.sin(wd*dt)
            C0=np.exp(-xi*w*dt)*np.cos(wd*dt)
            S1=-w_*S0+wd*C0
            C1=-w_*C0-wd*S0
            S2=-a1*S0-a2*C0
            C2=-a1*C0+a2*S0
            B=np.array([
            [S0, C0, 1, dt,dt**2,  dt**3],
            [S1, C1, 0,  1, 2*dt,3*dt**2],
            [S2, C2, 0,  0,    2,   6*dt]
            ])
            C=np.array([
            [-wd, -w_,   0,   1,     0,     0],
            [  0,   1,   1,   0,     0,     0],
            [  0,   0,w**2,  a0,     2,     0],
            [  0,   0,   0,w**2,  2*a0,     6],
            [  0,   0,   0,   0,2*w**2,  6*a0],
            [  0,   0,   0,   0,     0,6*w**2]
            ])
            C=np.linalg.inv(C)
            A=np.dot(B,C)
            y_=[y0]
            for i in range(len(T)-1):
                a=[y_[-1][1],y_[-1][0],R0[i],R1[i],R2[i],R3[i]]
                R_=np.array(a)
                y_.append(np.dot(A,R_))
            un=[]
            vn=[]
            an=[]
            for yt_ in y_:
                un.append(np.dot(mode,yt_[0]))
                vn.append(np.dot(mode,yt_[1]))
                an.append(np.dot(mode,yt_[2]))
            u.append(np.array(un))
            v.append(np.array(vn))
            a.append(np.array(an))
        us=sum(u)
        vs=sum(u)
        aas=sum(u)
        pd.Series(us[2]).plot()
        df=pd.DataFrame({'t':T,'u':us})
        
        
        
    def NewmarkBeta(self,T,F,u0,v0,a0,beta=0.25,gamma=0.5):
        """
        beta,gamma: parameters.\n
        u0,v0,a0: initial state.\n
        T: time list with uniform interval.\n
        F: list of time-dependent force vectors.
        """
        dt=T[1]-T[0]
        b1=1/(beta*dt*dt)
        b2=-1/(beta*dt)
        b3=1/2/beta-1
        b4=gamma*dt*b1
        b5=1+gamma*dt*b2
        b6=dt*(1+gamma*b3-gamma)
        K_=self.K+b4*self.C+b1*self.M 
        u=[u0]
        v=[v0]
        a=[a0]
        tt=[0]
        for (t,ft) in zip(T,F):
            ft_=ft+np.dot(self.M,b1*u[-1]-b2*v[-1]-b3*a[-1])+np.dot(self.C,b4*u[-1]-b5*v[-1]-b6*a[-1])
            ut=np.linalg.solve(K_,ft_)
            vt=b4*(ut-u[-1])+b5*v[-1]+b6*a[-1]
            at=b1*(ut-u[-1])+b2*v[-1]+b3*a[-1]
            u.append(ut)
            v.append(vt)
            a.append(at)
            tt.append(t)
        df=pd.DataFrame({'t':tt,'u':u,'v':v,'a':a})
        return df
        
    def WilsonTheta(self,T,F,u0=0,v0=0,a0=0,beta=0.25,gamma=0.5,theta=1.4):
        """
        beta,gamma,theta: parameters.\n
        u0,v0,a0: initial state.\n
        T: time list with uniform interval.\n
        F: list of time-dependent force vectors.
        """
        dt=T[1]-T[0]
        dt_=theta*dt
        b1=1/(beta*dt_**2)
        b2=-1/(beta*dt_)
        b3=(1/2-beta)/beta
        b4=gamma*dt_*b1
        b5=1+gamma*dt_*b2
        b6=dt_*(1+gamma*b3-gamma)
        K_=self.K+b4*self.C+b1*self.M 
        u=[u0]
        v=[v0]
        a=[a0]
        tt=[0]

        R_=[F[0]]
        for i in range(len(F)-1):
            R_.append(F[i]+theta*(F[i+1]-F[i]))
            
        for (t,ft) in zip(T,R_):
            ft_=ft+np.dot(self.M,b1*u[-1]-b2*v[-1]-b3*a[-1])+np.dot(self.C,b4*u[-1]-b5*v[-1]-b6*a[-1])
            ut_=np.linalg.solve(K_,ft_)
            vt_=b4*(ut_-u[-1])+b5*v[-1]+b6*a[-1]
            at_=b1*(ut_-u[-1])+b2*v[-1]+b3*a[-1]          
            
            at=a[-1]+1/theta*(at_-a[-1])
            vt=v[-1]+((1-gamma)*a[-1]+gamma*at)*dt
            ut=u[-1]+v[-1]*dt+(1/2-beta)*a[-1]*dt**2+beta*at*dt**2

            u.append(ut)
            v.append(vt)
            a.append(at)
            tt.append(t)
        df=pd.DataFrame({'t':tt,'u':u,'v':v,'a':a})
        return df
            
def Test():        
    sys = SDOFSystem(17.508,40*175,0.0)
    #harmonical load
    time=[]
    p=[]
    precise=[]
    A=10*4.448
    n=10
    T=10
    for i in range(T*500+1):
        t=i/500
        time.append(t)
#        precise.append(0.0254*0.33*(np.cos(10*t)-np.cos(20*t)))
        precise.append(0.0254*(0.32*np.cos(10*t-0.26)-np.exp(-4*t)*(0.31*np.cos(19.6*t)+0.11*np.sin(19.6*t))))
#        p.append(np.sin(2*np.pi/T*t*n)*A)
        p.append(A*np.cos(10*t))
        
#        p.append(0)
#        
#        if i<=100:
#            p.append(1)
#        else:
#            p.append(0)
        
#    df=pd.DataFrame({'t':time,'p':p})
#    df['p'].plot()
    
#    u=[]
#    for i in range(T*100):
#        t_list=time[:(i+1)]
#        p_list=p[:(i+1)]
#        u.append(sys.Duhamel(t_list,p_list))
#    
#        
#    df=pd.DataFrame({'t':time[1:],'u':u,'precise':precise[1:]})
#    df.plot('t','u')
#    df.plot('t','precise')
    
    df=sys.InterpExcite(time,p)
    df.plot('t','u')
    df2=sys.AverageAcceleration(time,p)
    df2.plot('t','u')
    
def MTest():
#    #Test with the mode-decomposite-method result in Zhu's Structural Mechanics p.134
#    m=2e5
#    k=200e6
#    t=np.arange(0.01,10,0.01)
#    A=0.1    
#    y=100*A*np.sin(10*t)
#    F=np.array([1.5*m*y,m*y,1.5*m*y]).T   
#    K=np.array([
#    [2*k,  -k,  0],
#    [ -k, 2*k, -k],
#    [  0,  -k,  k]
#    ])
#    M=np.array([
#    [1.5*m, 0,     0],
#    [    0, m,     0],
#    [    0, 0, 1.5*m]
#    ])
#    C=np.zeros((3,3))
#    msys=MDOFSystem(M,K,C)
#    u0=v0=a0=np.zeros(3)
#    df1=msys.NewmarkBeta(t,F,u0,v0,a0,beta=0.25,gamma=0.5)
#    df2=msys.WilsonTheta(t,F,u0,v0,a0,beta=0.25,gamma=0.5,theta=1.42)
#    df3=msys.ModalDecomposition(t,F,u0,v0,a0,xi=0)
#    return df1,df2,df3    
    m=1
    k=100e3
    K=np.array([
    [2*k,  -k,  0],
    [ -k, 2*k, -k],
    [  0,  -k,  k]
    ])
    M=np.array([
    [m, 0, 0],
    [0, m, 0],
    [0, 0, m]
    ])
    DOF=3
    C=np.zeros((DOF,DOF))
    msys=MDOFSystem(M,K,C)
    w,f,T,mode=msys.EigenMode(DOF)
    M_=np.dot(np.dot(mode.T,M),mode)#generalized mass
    px=[]
    Vx=[]
    Xm=[]
    mx=np.diag(M)
    for i in range(len(mode)):
        px.append(-np.dot(mode[:,i].T,mx))
        Vx.append(px[-1]**2)
        Xm.append(Vx[-1]/3/m)
    spec=(t,f)
    a=np.interp(T[i],t,f)
    
def ETest():
    #Test with the mode-decomposite-method result in Zhu's Structural Mechanics p.134
    m=2e5
    k=200e6
    t=np.arange(0.001,10,0.001)
    A=0.1    
    y=100*A*np.sin(10*t)
    F=np.array([1.5*m*y,m*y,1.5*m*y]).T   
    K=np.array([
    [2*k,  -k,  0],
    [ -k, 2*k, -k],
    [  0,  -k,  k]
    ])
    M=np.array([
    [1.5*m, 0,     0],
    [    0, m,     0],
    [    0, 0, 1.5*m]
    ])
    print(np.rank(K))
    C=np.zeros((3,3))
    msys=MDOFSystem(M,K,C)
    df=msys.EigenMode(1)
    return df

#if __name__=='__main__':
#    df1,df2,df3=MTest()
#    sdf1=pd.DataFrame(index=df1.t)
#    sdf2=pd.DataFrame(index=df2.t)
#    sdf3=pd.DataFrame(index=df3.t)
#    for j in range(3):
#        u0=[]
#        u1=[]
#        for i in df1['u']:
#            u0.append(i[j])
#        for i in df2['u']:        
#            u1.append(i[j])
#        for i in df3['u']:        
#            u1.append(i[j])
#        sdf1[str(j)]=u0
#        sdf2[str(j)]=u1
#    u3=[]
#    for i in df3['u']:        
#        u3.append(i[j])
#    sdf3['u']=u3
#    sdf1.plot()
#    sdf2.plot()
#    sdf3.plot()
##    df=ETest()
if __name__=='__main__':
    MTest()

