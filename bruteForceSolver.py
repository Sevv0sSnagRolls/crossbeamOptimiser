# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:39:01 2020

@author: Doug

this module is dedicated to the generation of descriptor variables for making trusses

"""

import numpy as np
import trussObjects
import TrussAssembly
import trussAssessment
import time
import math
from operator import itemgetter

def bruteForce(criterion, ho, wo, w, L, wt, ht, FOS, E, S_uc, rho, P, P_lat):
    '''
    iterates through all combinations
    
    Saves truss obejcts that work
    
    '''
    
    #initialise timer and count for the time
    t2 = int(round(time.time() * 1000))
    counter = 0
    
    #solutions container
    trusses = []
    
    #iterate through how many sections to evalute for the geometry
    for n in sections:      
        iterables = list( generateGlobalStructure(n, ho, wo, w, L, wt, ht) )
        totalIterationsEstimate = len(sections)*len(iterables) #variable for the timer
        
        #dumb line, but was easier to build constructor using discritisation due to offset of pyramid
        n = 2*n
        
        #iterate through each geometry type
        for variables in iterables:
            t2 = printTimeEstimate(t2, counter, totalIterationsEstimate)
            
            #generate the nodes and the initial guesses for the member sizing
            nodeObjectHandles, memberDict = createNodeandMemberObjects(n, variables, w, L, wt, ht, FOS, E, S_uc, rho, P, P_lat, K=0.5)
            
            #generate the load and support matricies for the stiffness matrix solver
            loadCases = loadsAndSupports(P, P_lat, len(nodeObjectHandles))

            #generate the truss object
            trussObjectHandle, success = TrussAssembly.trussAssembler(nodeObjectHandles, memberDict, 'Pyramidal')
            
            #if element matrix was created successfully etc
            if success:
                #test the truss against the loads
                #use the assessment tool to evaluate all the load cases against the truss
                trussPass = trussAssessment.assessTruss(  trussObjectHandle, loadCases, criterion )
                
                if trussPass:
                    trussObjectHandle.calculateMass()
                    trusses.append( [trussObjectHandle.mass, trussObjectHandle] )

                else:
                    """Once you have the geometry/global shap defined, need to iterate through the rod sizes until
                    the critera are met"""
                    
                    """NON STANDARD BAD FUNCTION THAT SHOULD BE A GRADIENT DSCENT, AGAIN TIME ISSUES"""
                    D_list = list( updateDiameters( memberDict['D'], d_step ) )
                    for D in D_list:
                        trussObjectHandle.updateFailedMembers(D)             
                        trussPass = trussAssessment.assessTruss(  trussObjectHandle, loadCases, criterion )
                        if trussPass:
                            trussObjectHandle.calculateMass()
                            trusses.append( [trussObjectHandle.mass, trussObjectHandle] )
                            break;
            counter += 1
        
    return trusses


def updateDiameters(D, d_step):
    """
    not a pretty function, works backwards as the primary members have the most effect on stiffness
    so will help find a solution faster in terms of ordering the produced combinations

    Parameters
    ----------
    D : TYPE
        DESCRIPTION.
    d_step : TYPE
        DESCRIPTION.

    Yields
    ------
    D : TYPE
        DESCRIPTION.

    """
    d1, d2, d3, d4, d5 = D
    for d5_n in listGenerator(d5, d5 + k*d_step, d_step):
        for d4_n in listGenerator(d4, d4 + k*d_step, d_step):
            for d3_n in listGenerator(d3, d3 + k*d_step, d_step):
                for d2_n in listGenerator(d2, d2 + k*d_step, d_step):
                    for d1_n in listGenerator(d1, d1 + k*d_step, d_step):
                        D = [d1_n, d2_n, d3_n, d4_n, d5_n]
                        yield D

       
def loadsAndSupports(P, P_lat, n):    
    supports = np.zeros( (n,3), dtype=int )
    supports[0] = [1,1,1]
    supports[1] = [1,1,1]
    supports[2] = [1,1,1]
    
    #downwards bending - flying a hull
    loadCases = []
    loadsCase1 = np.zeros( (n,3), dtype=float )
    loadsCase1[-3]= [0,0,-P/3]
    loadsCase1[-2]= [0,0,-P/3]
    loadsCase1[-1]= [0,0,-P/3]
    loadCases.append( [loadsCase1, supports] )
    
    #upwards bending - slam back into water
    loadsCase2 = np.zeros( (n,3), dtype=float )
    loadsCase2[-3]= [0,0,P/3]
    loadsCase2[-2]= [0,0,P/3]
    loadsCase2[-1]= [0,0,P/3]
    loadCases.append( [loadsCase2, supports] )
    
    #downwards bending - with a tilt angle to check any weird particulars with truss shape
    loadsCase3 = np.zeros( (n,3), dtype=float )
    deg = 25
    loadsCase3[-3]= [P*np.cos( (90-deg)*np.pi/180)/3, 0, P*np.cos(deg*np.pi/180)/3]
    loadsCase3[-2]= [P*np.cos( (90-deg)*np.pi/180)/3, 0, P*np.cos(deg*np.pi/180)/3]
    loadsCase3[-1]= [P*np.cos( (90-deg)*np.pi/180)/3, 0, P*np.cos(deg*np.pi/180)/3]
    loadCases.append( [loadsCase3, supports] )
    
    loadsCase4 = np.zeros( (n,3), dtype=float )
    loadsCase4[-3]= [0,0,-P_lat/3]
    loadsCase4[-2]= [0,0,-P_lat/3]
    loadsCase4[-1]= [0,0,-P_lat/3]
    loadCases.append( [loadsCase4, supports] )
        
    return loadCases



def printTimeEstimate(t2, counter, total):
    t1 = t2
    t2 = int(round(time.time() * 1000))
    dt = (t2 - t1) / 1000
    timeRemaining = round( dt * ( (total - counter) ) / 60, 2 )
    print( round( (counter)/(total)*100, 2), '% | remaining: ', timeRemaining, 'Minutes')
    return t2


def generateMembers(memberDict, d_step):
    dpp = [d + (d_step*(10**-3)) for d in memberDict['D'] ]
    diameters = dpp
    areas = [ np.pi*( (d/2)**2 ) for d in diameters]
    inertias = [ np.pi/4*( (d/2)**4 ) for d in diameters]
    memberDictNew = { 'D': diameters, 'A':areas, 'I':inertias, 'E':memberDict['E'], 'S_uc':memberDict['S_uc'], 'rho':memberDict['rho'] }
    return memberDictNew


def listGenerator(minVal, maxVal, stepSize):
    return np.linspace(minVal, maxVal, int( (maxVal - minVal)/stepSize) + 1)


def generateGlobalStructure(n, ho, wo, w, L, wt, ht):
    '''
    Designed to generate nodes in accordance with structure type
    
    A is the fixed bottom point of the mast
    
    B is the variable top of the truss at the mast
    
    C is the first point on the bottom of the truss at the height of the tender
    
    D is the second point of the bottom of the truss at the height of the tender
    
    E is the point above C at the top of the truss and is variable in height
    
    G is the point above D top of the truss and is variable in height
    
    F is the bottom of the truss at the connection to the outrigger
    
    H is the top of truss at connection to outrigger
    
    '''
    
    #this is the grid size hozirontally determined by the number of structural pieces
    s = L/(2*n)
    
    #first truss point
    B_z_StepSize = stepSizeVertical #[m]
    B_z_Min = 1.2 + h_h #[m] - fixed by Rob in email
    B_z_Max = 1.2 + h_h #[m] - fixed by Rob in email
    
    #disretisation is dumb, but
    C_x_StepSize = 2*s #[m]
    #got to start out atleast one unit away from the mast - no check here to see whether it's outside of hull width
    C_x_Min = s #[m]
    #how many points is the tender in width
    nst = math.ceil( wt/s )
    C_x_Max = (n - nst)*s
   
    E_z_StepSize = stepSizeVertical #[m]

    F_z_StepSize = stepSizeVertical #[m]
    F_z_Min = 0.4 #[m] - fixed by Rob in email
    F_z_Max = 0.4 #[m] - fixed by Rob in email
    
    #specfic generator for the points - not labelled very well - see diagram in my notes
    H_z_StepSize = stepSizeVertical #[m]    
    for B_z in listGenerator(B_z_Min, B_z_Max, B_z_StepSize):
        for C_x in listGenerator(C_x_Min, C_x_Max, C_x_StepSize):
            E_x = C_x 
            E_z_Min = B_z
            E_z_Max = B_z + maxTrussDeltaHeight  #[m]
            D_x = C_x + (nst + 1)*s
            for E_z in listGenerator(E_z_Min, E_z_Max, E_z_StepSize):
                G_x = D_x
                G_z = E_z
                for F_z in listGenerator(F_z_Min, F_z_Max, F_z_StepSize):
                    H_z_Min = F_z + minTrussDepth #[m]
                    H_z_Max = B_z_Max + H_z_StepSize #[m]  
                    for H_z in listGenerator(H_z_Min, H_z_Max, H_z_StepSize):                      
                        yield B_z, C_x, D_x, E_x, E_z, F_z, G_x, G_z, H_z    



def createNodeandMemberObjects(n, globalStructure, w, L, wt, ht, FOS, E, S_uc, rho, P, P_lat, K=0.5):
    '''
    Parameters
    ----------
    n : TYPE
        DESCRIPTION.
    globalStructure : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    wt : TYPE
        DESCRIPTION.
    ht : TYPE
        DESCRIPTION.
    FOS : TYPE
        DESCRIPTION.
    E : TYPE
        DESCRIPTION.
    S_uc : TYPE
        DESCRIPTION.
    rho : TYPE
        DESCRIPTION.
    P : TYPE
        DESCRIPTION.
    P_lat : TYPE
        DESCRIPTION.
    K : TYPE, optional
        DESCRIPTION. The default is 0.5.

    Returns
    -------
    nodes : TYPE
        DESCRIPTION.
    memberDict : TYPE
        DESCRIPTION.

    '''
    #using this vector, we can make the node objects
    B_z, C_x, D_x, E_x, E_z, F_z, G_x, G_z, H_z = globalStructure
    
    #using 2*n to make sure there is space for the offsets for the truss
    x_divisions = np.linspace(0, L, n)

    #iterate along x and generate the points
    i = 0
    nodes = []
    j = 0
    for x in x_divisions:
        
        if j == 0 or (j % 2) == 1 or x == x_divisions[-1]:
            if x < C_x:
                m_1 = (ht - 0 )/(C_x - 0)
                z_1 = m_1*x + 0
            elif x >= C_x and x <= D_x:         
                z_1 = ht
            elif x > D_x:
                m_1 = (F_z - ht)/(L - D_x)
                z_1 = m_1*(x - D_x) + ht
            nodes.append( trussObjects.node(i, np.array([x, w/2, z_1]) ) )
            i += 1
            
        if j == 0 or (j % 2) == 0 or x == x_divisions[-1]:
            if x < E_x:
                m_2 = (E_z - B_z)/(E_x - 0)
                z_2 = m_2*x + B_z
            elif x >= E_x and x <= G_x:
                z_2 = E_z
            elif x > G_x:
                m_2 = (H_z - G_z )/(L - G_x)
                z_2 = m_2*(x - G_x) + G_z
            
            nodes.append( trussObjects.node(i, np.array([x, 0, z_2]) ) )
            i += 1
            
            nodes.append( trussObjects.node(i, np.array([x, w, z_2]) ) )
            i += 1                                           
    
        j += 1
    
    
    x_StepSize = x_divisions[2] - x_divisions[1]    
    #need to size the elements here or give atleast some indication of a starting diameter
    l1 = 2*x_StepSize
    l2 = 2*x_StepSize
    l3 = ( (x_StepSize)**2 + (w/2)**2 + (B_z)**2 )**(0.5)
    l4 = ( (2*x_StepSize)**2 + (w)**2 )**(0.5)
    l5 = w
    
    #member angles
    theta3 = np.arccos(B_z/l3)
    theta4 = np.arccos(2*x_StepSize/l4)
    
    #member approximate axial loads
    F1 = P*L/( B_z )*FOS
    F2 = F1/2
    F3 = P / ( 2*np.cos(theta3) )*FOS
    F4 = P_lat/(np.cos(theta4))*FOS
    F5 = P_lat*FOS
    
    #member section sizing initialisation - the lazy way
    d1_S_uc = 2*( F1/S_uc*(1/np.pi) )**0.5
    d2_S_uc = 2*( F2/S_uc*(1/np.pi) )**0.5
    d3_S_uc = 2*( F3/S_uc*(1/np.pi) )**0.5
    d4_S_uc = 2*( F4/S_uc*(1/np.pi) )**0.5
    d5_S_uc = 2*( F5/S_uc*(1/np.pi) )**0.5
    
    d1_Buckle = 2*( F1*((l1*K)**2) / (0.85*E*(np.pi**3)/4) )**(1/4)
    d2_Buckle = 2*( F2*((l2*K)**2) / (0.85*E*(np.pi**3)/4) )**(1/4)
    d3_Buckle = 2*( F3*((l3*K)**2) / (0.85*E*(np.pi**3)/4) )**(1/4)
    d4_Buckle = 2*( F4*((l4*K)**2) / (0.85*E*(np.pi**3)/4) )**(1/4)
    d5_Buckle = 2*( F5*((l5*K)**2) / (0.85*E*(np.pi**3)/4) )**(1/4)
    
    d1 = max(d1_S_uc, d1_Buckle) - k2*d_step
    d2 = max(d2_S_uc, d2_Buckle) - k2*d_step
    d3 = max(d3_S_uc, d3_Buckle) - k2*d_step
    d4 = max(d4_S_uc, d4_Buckle) - k2*d_step
    d5 = max(d5_S_uc, d5_Buckle) - k2*d_step
    
    diameters  = [d1, d2, d3, d4, d5]  
    areas = [ np.pi*( (d/2)**2 ) for d in diameters]
    inertias = [ np.pi/4*( (d/2)**4 ) for d in diameters]
    E = [E, E, E, E, E]
    S_uc = [S_uc, S_uc, S_uc, S_uc, S_uc ]
    rho = [rho,rho,rho,rho,rho]
    
    memberDict = { 'D': diameters, 'A':areas, 'I':inertias, 'E':E, 'S_uc':S_uc, 'rho':rho }
    
    return  nodes, memberDict


'''all heights are measured from the waterline'''

#-----------------------------------------------------------------------------
#PROA DIMENSIONS
#Centre to Centre Length 
L_BCB = 8.2 #[m]

h_h = 0.3 #[m] - This is the freeboard of the hull full laden which will add to the height of the truss from the water

#Truss width
w = 0.5 #[m] the desired truss width


#tender dimesions 
ht = 1.0 #[m] this is the clearance desired between water surface and bottom of the truss
wt = 3.5 #[m] tender width plus a working margin

#outrigger depth/2
ho = 2.0/2
#outrigger width/2
wo = 1.2/2

radius_outrigger_stub_mast = 0.15 #[m]
radius_hull_stub_mast = 0.15 #[m]

#get the correct distance so then the truss entities can have the members go down to the correct load point.
# L = L_BCB - radius_outrigger_stub_mast - radius_hull_stub_mast
L = 9.0 #[m]- Set by Rob


m_hull = 13.0 #[ton] - full laden with structural mass estimate
m_outrigger = 2.0 #[ton] - max laden load in outrigger when flying a hull
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#LOADS
#Righting Moment
k_crossbeam_Stiffness_Path_Factor = 0.9
k_dynm = 3.0
P = m_outrigger*k_dynm*k_crossbeam_Stiffness_Path_Factor*9.81*(10**3) #[N]
#k_crossbeam_Stiffness_Path_Factor
#remember the tricky thing here is I'm loading all the outrigger mass through one crossbeam as
#there is the potential for only one mast to causea righting moment (one sail condition)
#and the hull has very little torisonal stiffness/ is a less stiff load path so one crossbeam will take far more load
#this is not an isolated issue though and changing the hull design will change the load path
#somthing to check later, but this is the safest assumption. 
#If it's too much mass, you can use the factor provided to alter the percentage load

#longitudinal/frontal impact load
k_F_LT_h = 2.5 #- 2.5g's for main hull deceleration from 12215-7 global loads
k_F_LT_o = 5.0 #- 5.0g's for outer hull deceleration F_LT from 12215-7 global loads
P_lat = max( (1/2)*m_hull*k_F_LT_h*9.81*10**3, (1/2)*m_outrigger*k_F_LT_o*9.81*10**3) #[N]
#with this load, it's actually shared between the two crossbeams so it's a half unlike the first load case
#may be slightly more on the crossbeam closer to where the impact was
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#Material Properties
E = 45*10**9 #[Pa] - From Spreadsheet of rod mechanical properties estimate. Lines up with experiments
S_uc = 350*10**6 #[Pa] - An accetable value for this failure mode, may wish to reduce, very unlikely
rho = 1800 #[kg/m^3]
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#SOLVER FACTORS


#This set the minimum distance allowable between the top and bottom - more of a manufacturing constraint
minTrussDepth = 0.25
#this is the discritisation verticall
stepSizeVertical = 0.05
#This is the max height the truss can rise to above the B_z height/height of the top truss at the mast
maxTrussDeltaHeight = 0.5 #[m]
    
#size of step for the diameters of the rod when iterating to find a structure that passes criteria
d_step = 0.5*(10**-3) #size of the step [m] of the truss
    
#this is the number of sections to evaluate/number of segments in the truss
sections = range(6, 7, 1)

#criteria, generally 12215 uses a FOS of 2 on everything
FOS = 2.0
criterion = {   
                'maxDisplacement': 0.02*L, 
                'FOS' : FOS
            }
k = 10  #this is the number of iterations of the d_step done for each rod when trying to find suitable members for the global geometry to resolve criteria
k2 = 2 #this is a reduction factor for the intial diameter guess on buckling, reason is that load estimates for initialisation are not ideal, but are close

'''For more factors relating to global geo go to geo generate function
made more sense to keep them in there'''
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#run 
ans = bruteForce(criterion, ho, wo, w, L, wt, ht, FOS, E, S_uc, rho, P, P_lat)
#sort the trusses by mass ascending
ans = sorted(ans, key=itemgetter(0))
#-----------------------------------------------------------------------------
#create viewing class
class solutionViewer():
    
    def __init__(self, ans):
        self.solutions = ans

    def view(self, i):
        self.solutions[i][1].plot(L)
        print('')
        print('MASS: [kg]')
        print( round(self.solutions[i][1].mass, 0) )
        nodes = []
        print('')
        print('NODE POSITIONS: [m]')
        for node in self.solutions[i][1].nodes:
            nodes.append(node.position)
            print( np.round(node.position, 2) )
        print('')
        print('DIAMETERS: [mm]')
        D = [round(1000*d, 1) for d in self.solutions[i][1].D ]
        print(D)

#-----------------------------------------------------------------------------
#Solution viewing
sol = solutionViewer(ans)

#just muck around in console with this command to view others
i = 0
sol.view(i)

#-----------------------------------------------------------------------------


















