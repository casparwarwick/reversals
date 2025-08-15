#*******************************************************************************
#Reversing the reversal
#*******************************************************************************
#Python parts for reversal utility 
#*******************************************************************************

#-------------------------------------
#1. Set-up
#-------------------------------------

import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE" # Needed to suppress conflicts. 
import numpy as np
from sfi import Data
from sfi import Macro
from sfi import Matrix

from scipy.optimize import minimize
from scipy.optimize import LinearConstraint

alpha_value=float(Macro.getLocal('alpha'))
use_theil = Macro.getLocal('theil') != ''

#-------------------------------------
#1.1 Define cost function
#-------------------------------------

if use_theil:
    # Use normalized Theil index
    def cost(l):

        #-------------------------------------
        # 1.1.1 Compute differences and normalize
        #-------------------------------------
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
        
        # Normalize differences to sum to 1
        maxdl = maxl - minl
        dl_normalized = [d/maxdl for d in dl]
        
        #-------------------------------------
        # 1.1.2 Compute Theil T index
        #-------------------------------------
        n = len(dl_normalized)
        mean = 1.0 / n  # Since normalized differences sum to 1
        
        theil = 0.0
        for val in dl_normalized:
            if val > 0:
                theil += (val / mean) * np.log(val / mean)
        theil = theil / n
        
        # Normalize by maximum possible Theil index ln(n)
        max_theil = np.log(n)
        normalized_theil = theil / max_theil
        
        exponent = 1/alpha_value
        cost = (normalized_theil)**exponent
        return cost
        
else:
    # Use general alpha-based cost function
    def cost(l):
        
        #-------------------------------------
        #1.1.1 Basics
        #-------------------------------------
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        #-------------------------------------
        #1.1.2 dl
        #-------------------------------------
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
        maxdl = maxl - minl
        N = K - 1
        maxvar = (1/N - 1/N**2)*maxdl**2
        Edl = maxdl/N

        #-------------------------------------
        #1.1.3 Compute each component
        #-------------------------------------
        var_comp = [0]*(N)
        for i in range(0,N):
            var_comp[i] = (dl[i] - Edl)**2 
    
        #-------------------------------------
        #1.1.4 Get the variance, normalise, and output
        #-------------------------------------
        var = (1/N)*np.sum(var_comp)
        
        # Apply alpha parameter: cost = (var/maxvar)^(1/alpha)
        exponent = 1/alpha_value
        cost = (var/maxvar)**exponent
            
        return cost

#-------------------------------------
#1.2 Import coefficients, number of labels, scale_min and scale_max, and original labels from Stata
#-------------------------------------

#Coefficients
bd = Data.get(var=["bd"])
sign_var = Data.get(var=["sign"])
sign = sign_var[0]

#Reversal point
reversal_point_var = Data.get(var=["revpoint"])
reversal_point = reversal_point_var[0]

#Number of labels
nlabs = len(bd)  

#scale_min and scale_max
scale_min = float(Macro.getLocal('scale_min'))
scale_max = float(Macro.getLocal('scale_max'))

#Original labels
labels_depvar = np.asarray(Matrix.get("_labels_depvar"))[0:]
labels_depvar = labels_depvar.tolist()
l_original    = [val[0] for val in labels_depvar]
l_transformed = [val[0] for val in labels_depvar] # for an initial value

#-------------------------------------
#1.3 Set constraint that labels are weakly increasing
#-------------------------------------
	
#Define monotonicity constraints
#This will create a matrix with K-1 rows and K columns, with each row having 1 at position i and -1 at position i+1
monotone_array1 = []
for i in range(0,nlabs-1):
	j = i+1
	y = [0]*nlabs
	y[i] = 1
	y[j] = -1	
	monotone_array1.append(y) 
monotone_array2 = [-np.inf]*(nlabs-1)
monotone_array3 = [0]*(nlabs-1)

monotonicity_constraint = LinearConstraint(monotone_array1, monotone_array2, monotone_array3)

#-------------------------------------
#1.4 Set constraint that new labels need to lead to a reversal
#-------------------------------------

reversal_array = [0]*nlabs
for i in range(0,nlabs):
	h = i - 1
	if i==0:
		reversal_array[i] = bd[i]	
	if i>0 & i<(nlabs-1):
		reversal_array[i] = bd[i]-bd[h]
	if i==(nlabs-1):
		reversal_array[i] = -bd[h]

# Check if coefficient should cross the reversal_point from above or below
# If original sign is positive, we want the transformed coefficient to be <= reversal_point
# If original sign is negative, we want the transformed coefficient to be >= reversal_point		
if sign > 0:
	reversal_constraint = LinearConstraint(reversal_array, [-np.inf], [reversal_point])
else:
	reversal_constraint = LinearConstraint(reversal_array, [reversal_point], [np.inf])
	
#-------------------------------------
#1.5 Set constraint that labels need to be between 1 and the width of the scale, given by "width". 
#-------------------------------------

tmp1 = [0]*nlabs
tmp1[0]=1
tmp2 = [0]*nlabs
tmp2[nlabs-1]=1

boundary_constraint = LinearConstraint([tmp1,tmp2], [scale_min,scale_max], [scale_min,scale_max])

#-------------------------------------
#1.6 Minimize cost function subject to the constraints
#-------------------------------------

#Minimize cost function and save result
result = minimize(cost, l_transformed, constraints=[monotonicity_constraint, reversal_constraint, boundary_constraint])

#Save cost value
cost_value = [result.fun]*nlabs # just puts things into the right format 			

#-------------------------------------
#1.7 Output result to Stata
#-------------------------------------

Data.addVarDouble("python_labels")
Data.addVarDouble("python_cost")

Data.store("python_labels", None, result.x, None)
Data.store("python_cost", None, cost_value, None)