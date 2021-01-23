# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 10:27:12 2020

@author: Doug


This module is dedicated to generating truss objects 
and the helper functions that assist with building up truss objects


TO DO:
- Need to make trussTypeDict a little bit more wideused accross functions etc

"""


import numpy as np
import trussObjects


def trussAssembler(nodeObjectHandles, memberDicts, trussType, elementJoins=None):
    '''
    This is designed to generate truss objects 
    
    It accepets the base variables for the truss configuration:
    - Nodes and node positions
    - Members (variable, constant)
    - Truss connection style (user defined/standard)
        - standard joins can be generated by the elementJoinSequqenceGenerator()
    
    It will then put these all together into a truss object which can passes to
    a feasbility and/or evaluation function
    
    When it passes, the optimiser algorithm can save the truss object etc
    
    Basically it's designed to have a vector of variable inputs from any
    sort of optimiser and will build up the truss object
    
    This is one way to start a truss, then 
    '''
    
    n = len(nodeObjectHandles)
    
    #needs to be better optimised and utliused across functions
    trussTypeDict = { 
                      'Pyramidal': {'nP':3, 'nM':5}, 
                      'Cubic':     {'nP':4, 'nM':3}, 
                    }
    #!!!have not correctly checked the cubic setup - might need lots of debugging!!
    
    #correctly setup for the element join generator
#    if trussType != 'custom':
#        if trussTypeDict[trussType]['nM'] < len(memberDicts):
#            variableMembers = True
#        else:
#            variableMembers = False
    
    elementJoins = elementJoiningSequenceGenerator(n, trussType, variableMembers=False)

    elementObjectHandles, success = createElementObjects(elementJoins, nodeObjectHandles, memberDicts)
    
    if success:
        #create the truss object
        trussObjectHandle = trussObjects.truss(nodeObjectHandles, elementObjectHandles, memberDicts['D'])    
        return trussObjectHandle, success
    
    return None, success



def createElementObjects(elementJoins, nodeObjectHandles, memberDicts):
    '''
    Runs on the assumption that nodes are integer labels and it's done correctly etc
    Also runs on the assumption that the members ID matrix lines up correctly
    members dict needs to be correctly ordered and have the correct properties to create
    the elements
    
    Input 
    - joins matrix
    - members Dict
    
    Outputs
    - Element Objects
    - Whether operation was successful
    
    
    '''
    elementObjectHandles = []
    successFlag = True
    n, m = np.shape(elementJoins)
    for i in range(1, n):
        for j in range(0, i):
            value = elementJoins[i][j]
            if value > 0:
                nodei = nodeObjectHandles[i]
                nodej = nodeObjectHandles[j]
                #generate an element object and add handle to list of handles
                elementObjectHandles.append( trussObjects.element(
                                             nodei, 
                                             nodej, 
                                             value,
                                             memberDicts['A'][value-1], 
                                             memberDicts['E'][value-1], 
                                             memberDicts['I'][value-1], 
                                             memberDicts['S_uc'][value-1], 
                                             memberDicts['rho'][value-1], 
                                             memberDicts['D'][value-1]) 
                                            )
                    
    return elementObjectHandles, successFlag




def elementJoiningSequenceGenerator(n, trussType, variableMembers=False):
    '''
    This function is for generating standard types of trusses to evaluate
    The user can choose to generate more random connection styles in their own function
    and by pass this 'helper function if they choose to
    The default accepted stucture of the n x n node connection matrix is important
    
    Inputs
    n = number of nodes in the truss design
    trussType = [pyramidal, square...]
    
    Function will then generate the element joining array
    for the specific truss type of size n x n
    
    variable members will allow member section attrbutes ID's to increment 
    with each section generated. This means that a 12 node, 3 section pyramidal truss
    can have different cross section members in each of the sections
    DETAIL
     - calculates the maximum number of expected truss members and chops the 
        end to make sure you don't have nM extras member section types 
    
    '''
    trussTypes = ['Pyramidal', 'Cubic']
    
    assert trussType in trussTypes
    assert type(n) is int
    
    #create the empty joins array of size n x n
    joins = np.zeros((n,n), dtype=int)
    
    #follow the specific generation sequence for the style of joining desired
    if trussType == 'Pyramidal':
        '''
        Pyramidal Truss Style
        
        Refer to docs if more is desired to be known about the sequencing
        '''
        nP = 3 #number of primary members -> This is the main driver of the repeated pattern size
        nM = 5 #number of types of members in the truss (which will mean the assignment of the joins and what joins them)
        
        #this is a caveat which allows each member to be different in the truss
        #so we can assign up to nP*nM different
        if variableMembers:
            maxVariableMemberValue = nP*nM
            
        pattern1 = np.array([ [3],
                              [3],
                              [1] ])
        
        pattern2 = np.array([ [5],
                              [3],
                              [2],
                              [4] ])
        
        pattern3 = np.array([ [3],
                              [0],
                              [2] ])
    
        patterns = [pattern1, pattern2, pattern3]
        
        for i in range(0, n-1):
            index = i % len(patterns)
            pat = patterns[index].copy()
            mpt, npt = np.shape(pat)
            
            #function for adding variable member indicies/ID's 
            #to the generation of the truss connection matrix
            if variableMembers:
                offset = nM*int(i/nP)
                if offset > (maxVariableMemberValue - nM):
                    offset -= nM
                pat += offset
                
            for j in range(0, mpt):
                #try catch for a very lazy method of overlapped matricies
                try:
                    joins[i+1 + j][i] = pat[j][0]
                except:
                    IndexError
        
        #key line to make sure the end connects together
        #dependent on the node generator working well
        if n == 14:
            joins[n-1][n-3] = 3
            
    return joins






if __name__ == "__main__":    
    #TEST THE FUNCTIONS
    testEJSG_1 =  elementJoiningSequenceGenerator(14, 'Pyramidal', variableMembers=False)
    print(testEJSG_1)


                        



















