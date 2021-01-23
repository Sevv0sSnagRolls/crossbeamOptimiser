# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 15:07:49 2020

@author: Doug

Module is dedicated to assessment of the truss structure against all load cases
it needs to be designed for

"""

import numpy as np
import FEMSolver
import trussObjects
import TrussAssembly

def assessTruss( trussObjectHandle, loadCases, criterion ):
    '''
    returns True if passes
    returns False is fails
    
    LoadCases Structure
    [ [ loads Matrix, supports matrix ], [ loads Matrix, supports matrix ], ... ]
    
    Where each matrix is the size of the truss (number of nodes)*3 dof to restrain
    As only dealing with axial forces
    '''
    
    for loadCase in loadCases:
        loads = loadCase[0]
        
        supports = loadCase[1]
       
        U, F, success = FEMSolver.solveTruss(trussObjectHandle, loads, supports)       
        
        trussPass = False
        #If solver finds a possible solution - check all the conditions specified
        #reason solver would not find a solution include a signular matrix from
        #ill specified setup, improper supports or strange load applications.
        if success:
            trussObjectHandle.checkStructuralFailure(U, criterion['FOS'])
                
            #need to write a smarter routine to check truss failure.
            #inbuilt method of the truss class perhaps.
            if trussObjectHandle.structuralFailure:
                trussPass = False
                
            elif trussObjectHandle.checkMaxDisplacement(U, criterion['maxDisplacement']):
                trussPass = False
                
            else:
                trussPass = True
                
        else:
            success = False
            break;
    
    if trussPass:
        return True
    
    return False





















#
#if __name__ =="__main_":
#    
#    L = 8.2
#
#    criterion = {   
#                    'maxDisplacement': 0.05*L, 
#                    'FOS' : 2.0
#                }
#    












