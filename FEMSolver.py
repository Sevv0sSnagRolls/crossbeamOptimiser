# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:46:33 2020

@author: Doug

Module is dedicated to Global stuffness matrix axial loads 3D solver

"""

import numpy as np

def globalStiffnessMatrix(elements, n, dof):
    #global stiffness matrix
#    global Kp
    
    Kp = np.zeros((dof, dof), dtype=float)

    for element in elements:
        indexA = 3*element.nodeA_ID
        indexB = 3*element.nodeB_ID
        
        A = element.A
        E = element.E    
        
        L = element.L     
        Cx = element.Cx 
        Cy = element.Cy 
        Cz = element.Cz 
        
        
        #lambda matrix
        K_lambda = (A*E/L)*np.array( [ [ Cx**2, Cx*Cy, Cx*Cz ],
                                       [ Cx*Cy, Cy**2, Cy*Cz ],
                                       [ Cx*Cz, Cy*Cz, Cz**2 ] ] )
        
        lx, ly = np.shape( K_lambda )
        
        #Lambda matrix AA global stiffness coordinates
        KAA_x1 = indexA
        KAA_x2 = indexA + lx
        KAA_y1 = indexA
        KAA_y2 = indexA + ly
        
        #Lambda matrix AB global stiffness coordinates
        KAB_x1 = indexA
        KAB_x2 = indexA + lx
        KAB_y1 = indexB
        KAB_y2 = indexB + ly
        
        #Lambda matrix BA global stiffness coordinates
        KBA_x1 = indexB
        KBA_x2 = indexB + lx
        KBA_y1 = indexA
        KBA_y2 = indexA + ly
        
        #Lambda matrix BB global stiffness coordinates
        KBB_x1 = indexB
        KBB_x2 = indexB + lx
        KBB_y1 = indexB
        KBB_y2 = indexB + ly
        
        #Add lambda matricies to global stiffness matrix
        #AA (+)
        Kp[KAA_x1:KAA_x2, KAA_y1:KAA_y2] += K_lambda
        #AB (-)
        Kp[KAB_x1:KAB_x2, KAB_y1:KAB_y2] -= K_lambda
        #BA (-)
        Kp[KBA_x1:KBA_x2, KBA_y1:KBA_y2] -= K_lambda
        #BB (+)
        Kp[KBB_x1:KBB_x2, KBB_y1:KBB_y2] += K_lambda
    
    return Kp
 
    
def solveDisplacements(nodes, elements, loads, supports, dof, Kp):
    #STEP 1 - Solve for displacements
    
    #set the force vector to the applied forces
    F = np.transpose( loads )
    
    #set the rows and columns of the stiffness matrix = 0 where displacements are 0
#    global Kp_z 
    Kp_z = np.copy(Kp)
    for i in range(0, dof):
        if supports[i]:
            #set rows and columns equal to zero associated with the node that is restrained
            Kp_z[i, :] = 0
            Kp_z[:, i] = 0
            
            #set the diagonal equal to 1
            Kp_z[i, i] = 1
            
            F[i] = 0
            
    #solve for nodal displacements
    #uses pinv which can have errors with singular stiffness matricies etc
    #try and catch the expection and return a failed result for the 
    #user to change the input params
    try:
        U = np.dot( np.linalg.pinv(Kp_z), F )
        success = True
    except Exception:
        U = 0
        success = False
    
    return U, success


def solveTruss(truss, loads, supports):   
       
    nodes = truss.nodes
    
    elements = truss.elements
    
    n = len(nodes)
    
    dof = 3*n
    
    nodes3D = []
    for node in nodes:
        nodes3D.append( node.position[0] )
        nodes3D.append( node.position[1] )
        nodes3D.append( node.position[2] )
    
    loads3D = []
    for load in loads:
        loads3D.append( load[0] )
        loads3D.append( load[1] )
        loads3D.append( load[2] )
        
    supports3D = []
    for support in supports:
        supports3D.append( support[0] )
        supports3D.append( support[1] )
        supports3D.append( support[2] )
    
    Kp = globalStiffnessMatrix(elements, n, dof)

    #STEP 1 - Solve for Displacements
    U, success = solveDisplacements(nodes3D, elements, loads3D, supports3D, dof, Kp)

    #STEP 2 - Solve for Reactions
    if success == True:
        F = np.dot(Kp, U)    
    else:
        F = 0

    return U, F, success





















