# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 14:10:53 2020

@author: Doug

This module is dedicated to the truss object classes

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


class truss():    
    
    """
    Truss Class
    
    Information
    Hold all the information of a generated truss object under the world conditions given
    this includes:
    
    Point A, Point B
    
    Nodes
    - Node Object Handles
    
    Elements
    - Element Connection Matrix (How to assemble the truss)
    - Element Objects Handles List
    
    Geometric Contraints
    
    Loadings
    
    
    
    .....
    
    
    Methods:
    Plot()
    Allows visualisation of the truss object
    
    calculateMass()
        
    

    """
    
    def __init__(self, nodes, elements, diameters):
        self.nodes = nodes
        self.elements = elements
        self.n = len(self.nodes)
        self.D = diameters
        
    def plot(self, L):
        '''
        SHould write a method to get min and max truss x node positions to
        get the length to callibrate the field of nodes to set aspect ratio correct
        ;'''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for element in self.elements:
            x = [element.xA, element.xB]
            y = [element.yA, element.yB]
            z = [element.zA, element.zB]
            ax.plot(x, y, z)
        max_range = np.array([L, L/4, L/4]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(L)
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(L/4)
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(L/4)
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w') 
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('z [m]')
        ax.set_title('Truss')
        return
        
    def calculateMass(self):
        self.mass = np.sum( [ element.findMass() for element in self.elements ] )
        return
    
    def checkStructuralFailure(self, U, FOS):
        '''
        runs the element calculation for failure on each element
        
        tells you if it fails under the current loads etc
        '''
        self.structuralFailure = False
        for element in self.elements:
            if element.checkStructuralFailure(U, FOS):
                self.structuralFailure = True
                return True
        return False
    
    def checkMaxDisplacement(self, U, critieria):
        '''
        runs the element calculation for failure on each element
        
        tells you if it fails under the current loads etc
        '''
        #group points back into (x,y,z)
#        U_array = np.array( [ np.linalg.norm( [ U[i], U[i+1], U[i+2] ] ) for i in range(0, len(U), 3) ] )
        self.maxDisplacement = np.max( abs(U) )
        if self.maxDisplacement <= critieria:
            return False
        return True
        
    def updateFailedMembers(self, D):
        for element in self.elements:
            element.updateDiameter(D[element.memberType - 1])
        self.D = D
        return 
        
    
class node():
    '''
    node object
    
    Setup
    Has an integer ID value for stiffness matrix assembly
    Has a numpy vector for position of node0 global 3D coordinates
    
    Methods
    Forces
    
    
    '''    
    def __init__(self, ID, position):
        self.ID = ID
        self.position = position
        
    def forces(self, F):
        self.force = np.array(F[0], F[1], F[2])


class element():
    '''
    node object
    
    Has an integer ID value for stiffness matrix assembly
    
    Has a numpy vector for position of node0 global 3D coordinates
    '''    
    def __init__(self, nodeA, nodeB, value, A, E, I, S_uc, rho, D):
        self.nodeA_ID = nodeA.ID
        self.nodeB_ID = nodeB.ID
        
        self.A = A
        self.I = I
        self.E = E
        self.S_uc = S_uc
        self.rho = rho
        
        self.D = D
        
        self.memberType = value
        
        positionA = nodeA.position
        positionB = nodeB.position
        self.xA = positionA[0]
        self.yA = positionA[1]
        self.zA = positionA[2]
        
        self.xB = positionB[0]
        self.yB = positionB[1]
        self.zB = positionB[2]
        
        self.dx = self.xB - self.xA
        self.dy = self.yB - self.yA
        self.dz = self.zB - self.zA
        
        self.L = ( (self.dx)**2 + 
                   (self.dy)**2 + 
                   (self.dz)**2   )**0.5
                  
        self.Cx = self.dx/self.L
        self.Cy = self.dy/self.L
        self.Cz = self.dz/self.L
        
        
    def internalForces(self, U):
        '''
        Solve for the internals based off of the solved displacement vector
        
        This could definitely be optimised from the truss parent object
        '''
        self.u = np.hstack( (U[3*self.nodeA_ID : (3*self.nodeA_ID + 3)], U[3*self.nodeB_ID : (3*self.nodeB_ID + 3)] ) ) 
        self.T = np.array([ [self.Cx, self.Cy, self.Cz, 0, 0, 0 ],
                            [ 0, 0, 0, self.Cx, self.Cy, self.Cz] ] )
        self.u_dash = np.matmul(self.T, self.u)
        self.internalForce = self.E*self.A*(self.u_dash[1] - self.u_dash[0])/self.L
        self.stress = self.internalForce/self.A
        return

    def checkStructuralFailure(self, U, FOS):
        #solve for the forces
        self.internalForces(U)
        #check the stress
        self.stressCheck(FOS)
        #check for buckling
        self.eulerBuckle(FOS)
        if self.compressiveFailure == True or self.bucklingFailure == True:
            self.fail = True
        else:
            self.fail = False
        
    def stressCheck(self, FOS):
        if self.stress < self.S_uc/FOS:
            self.compressiveFailure = False
        else:
            self.compressiveFailure = True
        return
            
    def eulerBuckle(self, FOS):
        K = 1.0 #assuming pin joints
        self.P_crit = 0.85*self.E*(np.pi**2)*self.I/( (self.L*K)**2 )
        #get abs to get if it's in compression
        if self.P_crit > abs( FOS*self.internalForce ):
            self.bucklingFailure = False
        else:
            self.bucklingFailure = True
        return
    
    def findMass(self):
        self.mass = self.L*self.rho*self.A
        return self.mass
    
    def updateDiameter(self, dpp):
        self.D += dpp
        self.A = np.pi*(self.D/2)**2
        self.I = np.pi/4*(self.D/2)**4
        if self.memberType == 1:
            self.A = 2*self.A
            self.I = 2*self.I











