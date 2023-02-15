import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
import nglview as nv
from nglview.contrib.movie import MovieMaker
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import dist
import matplotlib.ticker as tkr
import matplotlib.pylab as plb
import math as math
import MDAnalysis.transformations as trn
import scipy as sci
#this First calculates the the inertial moments of the protein, and then uses Singular Value Decomposition to calculate the eigenvectors corresponding xyz principal axes of the protein of interest.


# This decides what values to modify your plot xvalues by, if ommitted it will simply plot the index of the result array. 
def numfmt(x, pos): # your custom formatter function: divide by 10
    s = '{}'.format(x / 10)
    return s
yfmt = tkr.FuncFormatter(numfmt)


u = mda.Universe("sys.psf", "wrapprod.dcd")
import math as math
ag = u.select_atoms("(protein and name CA) and not resid 170:180 and not resname CYSF")
#This function calculates the principal axes of the protein by first calculating the inertial moments of the protein, and then Decomposing the resultant matrix to yield unit vectors corresponding to the xyz principal axes of the protein. For Simplicity, this only uses the alpha carbon atoms of each protein residue, can obviously be modified to use all protein backbone atoms
#It is very important to note that this calculation is technically for a rigid body, and therefore the more flexible your protein is, the less accurate the principal axis calculation will be. 
def calc_protein_principal_axis(ag):  
    com = ag.center_of_mass()
    mass = 12.011
    possubcom = ag.positions - com
    sqrmatrix = np.square(possubcom)
    ixx = -mass * (sqrmatrix[:,1].sum() + sqrmatrix[:,2].sum())
    ixy = mass * (possubcom[:,0] * possubcom[:,1]).sum()
    ixz = mass * (possubcom[:,0] * possubcom[:,2]).sum()
    iyy = -mass * (sqrmatrix[:,0].sum() + sqrmatrix[:,2].sum())
    iyz = mass * (possubcom[:,1] * possubcom[:,2]).sum()
    izz = -mass * (sqrmatrix[:,0].sum() + sqrmatrix[:,1].sum())
    inarray = np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])
    Eg, diag, prin = sci.linalg.svd(inarray)
    #principal_axis_array = np.array(resultarray[2])
    return prin
def vecmag(vector):
    return math.sqrt(sum(pow(element,2)for element in vector))
alphas = []
betas = []
gammas = []

#This part of the script then calculates the principal axes at each frame, and then the Tait-Bryan euler angles of the protein (Tilt, Spin, Rotation). Because of the way MDAnalysis positions the protein, The Principal Axis vectors must be corrected for direction, which the negative value check does. This was checked with a similar calculation done using the Orient plugin in VMD along with the Linear Algebra plugin. 
for ts in u.trajectory[::10]:  
    prinaxes = calc_protein_principal_axis(ag)
    xprinaxis = prinaxes[0]
    if xprinaxis[0] < 0:
        xprinaxis = xprinaxis * -1 
    yprinaxis = prinaxes[1]
    if yprinaxis[1] < 0:
        yprinaxis = yprinaxis * -1
    zprinaxis = prinaxes[2]
    if zprinaxis[2] < 0:
        zprinaxis = zprinaxis * -1
    print(zprinaxis)
    zpaxis_vector_length = np.linalg.norm(zprinaxis)
    a = (zprinaxis[2]) / zpaxis_vector_length
    #if a[2] < 0:
     
      #  zprinaxis = zprinaxis * -1
    b = np.array([0, 0, 70]) / np.linalg.norm(np.array([0, 0, 70]))
    beta = np.arccos(a)
    alpha = math.atan2(zprinaxis[0], -(zprinaxis[1]))
    alphas.append(alpha)
    gamma = math.atan2(xprinaxis[2], yprinaxis[2])
    gammas.append(gamma)
    betas.append(beta)
betaarray = np.degrees(np.array(betas))
alphaarray = np.degrees(np.array(alphas))
gammaarray = np.degrees(np.array(gammas))

