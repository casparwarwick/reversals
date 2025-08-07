#*******************************************************************************
#Reversing the reversal
#*******************************************************************************
#Python routine for MRS (coefficient ratio) reversal analysis
#*******************************************************************************

#=====================================
#1. Set-up
#=====================================

#-------------------------------------
#1.1 Import libraries
#-------------------------------------

import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

import numpy as np
from sfi import Data, Macro, Matrix
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, BFGS

#=====================================
#2. Define cost function (same as other scripts)
#=====================================

def cost(l):
    """
    Cost function for scale transformations
    Uses either alpha-parameterized or Theil index approach
    """
    cost_function_type = Macro.getLocal('cost_function_type')
    
    if cost_function_type == 'theil':
        #-------------------------------------
        #2.1 Theil index cost function
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
        
        # Compute Theil T index
        n = len(dl_normalized)
        mean = 1.0 / n  # Since normalized differences sum to 1
        
        theil = 0.0
        for i in range(n):
            if dl_normalized[i] > 0:
                theil += dl_normalized[i] * np.log(dl_normalized[i] / mean)
        
        # Normalize by maximum possible Theil index
        max_theil = np.log(n)
        cost_value = theil / max_theil
        return cost_value
        
    else:
        #-------------------------------------
        #2.2 Alpha-parameterized cost function (default)
        #-------------------------------------
        alpha_value = float(Macro.getLocal('alpha'))
        
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
            
        maxdl = maxl - minl
        N = K - 1
        Edl = maxdl/N

        var_comp = [0]*(N)
        for i in range(0,N):
            var_comp[i] = (dl[i] - Edl)**alpha_value 

        var = (1/N)*np.sum(var_comp)
        
        if alpha_value == 2:
            maxvar = (1/N - 1/N**2)*maxdl**alpha_value
        else:
            # For general alpha, compute maximum variance
            # Maximum occurs when all mass is at endpoints
            maxvar = 0.5 * (maxdl - Edl)**alpha_value + 0.5 * (0 - Edl)**alpha_value
        
        cost_value = (var/maxvar)**(1/alpha_value) if alpha_value != 1 else var/maxvar
        return cost_value

#=====================================
#3. Import data from Stata
#=====================================

#-------------------------------------
#3.1 Import coefficients and setup
#-------------------------------------

# Get denominator coefficients 
bdn = np.asarray(Matrix.get("_denominator_coeffs"))[0:]
bdn = bdn.tolist()
bdn = [val[0] for val in bdn]
bdn = np.asarray(bdn)

# Get numerator coefficients (matrix with each column being a different variable)
bdm_matrix = np.asarray(Matrix.get("_numerator_coeffs"))
num_vars = bdm_matrix.shape[1]

# Check if denominator is reversible
denom_reversible = int(Macro.getLocal('denom_reversible_flag')) == 1

# Number of labels
nlabs = len(bdn) + 1

# Original labels
l_original = np.array(range(1, nlabs+1, 1))

# Scale bounds
scale_min = np.amin(l_original)
scale_max = np.amax(l_original)

#=====================================
#4. Set up constraints
#=====================================

#-------------------------------------
#4.1 Monotonicity constraint
#-------------------------------------

monotone_array1 = []
x = [0]*nlabs
for i in range(0,nlabs-1):
    j = i+1
    y = x.copy()
    y[i] = 1
    y[j] = -1	
    monotone_array1.append(y) 
monotone_array1.append(x) 
monotone_array2 = [-np.inf]*nlabs
monotone_array3 = [0]*nlabs

monotonicity_constraint = LinearConstraint(monotone_array1, monotone_array2, monotone_array3)

#-------------------------------------
#4.2 Boundary constraint
#-------------------------------------

tmp1 = [0]*nlabs
tmp1[0]=1
tmp2 = [0]*nlabs
tmp2[nlabs-1]=1

boundary_constraint = LinearConstraint([tmp1,tmp2], [scale_min,scale_max], [scale_min,scale_max])

#=====================================
#5. Define ratio calculation functions
#=====================================

#-------------------------------------
#5.1 Function to calculate coefficient from transformation
#-------------------------------------

def calculate_coefficient(x, coeff_vector):
    """Calculate coefficient from transformation x and coefficient vector"""
    value = (x[0]-x[1])*coeff_vector[0]
    for i in range(1, nlabs-1):
        j = i+1
        value += (x[i]-x[j])*coeff_vector[i]
    return value

#-------------------------------------
#5.2 Function to calculate ratio
#-------------------------------------

def calculate_ratio(x, bdm_col):
    """Calculate ratio of numerator to denominator"""
    numer = calculate_coefficient(x, bdm_col)
    denom = calculate_coefficient(x, bdn)
    return numer/denom

#-------------------------------------
#5.3 Target ratio constraint function
#-------------------------------------

def ratio_constraint_func(x, bdm_col, target_ratio):
    """Constraint function for achieving target ratio"""
    current_ratio = calculate_ratio(x, bdm_col)
    return current_ratio - target_ratio

#=====================================
#6. Compute results for each numerator variable
#=====================================

#-------------------------------------
#6.1 Calculate original ratios and bounds
#-------------------------------------

# Store results for each variable
original_ratios = []
min_ratios = []
max_ratios = []
target_costs = []

for var_idx in range(num_vars):
    bdm_col = bdm_matrix[:, var_idx]
    
    #-------------------------------------
    #6.1.1 Calculate original ratio
    #-------------------------------------
    orig_ratio = calculate_ratio(l_original, bdm_col)
    original_ratios.append(orig_ratio)
    
    #-------------------------------------
    #6.1.2 Find min and max ratios from boundary transformations
    #-------------------------------------
    
    if denom_reversible:
        # If denominator is reversible, ratios are unbounded
        min_ratio = -np.inf
        max_ratio = np.inf
    else:
        # Get all possible ratios from d-regressions (boundary transformations)
        ratios = bdm_col / bdn
        min_ratio = np.amin(ratios)
        max_ratio = np.amax(ratios)
    
    min_ratios.append(min_ratio)
    max_ratios.append(max_ratio)
    
    #-------------------------------------
    #6.1.3 Calculate target cost (if target ratio specified)
    #-------------------------------------
    
    has_target = int(Macro.getLocal('has_target_ratio'))
    if has_target == 1:
        target_ratio = float(Macro.getLocal('target_ratio_value'))
        
        # Calculate target cost regardless of denominator reversibility
        # Note: if denominator is reversible, results may not be accurate but we still compute them
        
        # For reversible denominators, we still attempt the calculation 
        # but bounds checking is different (infinite bounds mean any ratio is theoretically achievable)
        if denom_reversible or (min_ratio <= target_ratio <= max_ratio):
            # Define constraint for target ratio
            target_constraint = NonlinearConstraint(
                lambda x: ratio_constraint_func(x, bdm_col, target_ratio), 
                0, 0, jac='2-point', hess=BFGS()
            )
            
            # Initial guess
            l_initial = np.random.uniform(low=l_original[0], high=l_original[nlabs-1], size=nlabs)
            l_initial = np.sort(l_initial)
            l_initial[0] = l_original[0]
            l_initial[nlabs-1] = l_original[nlabs-1]
            
            # Minimize cost subject to target ratio constraint
            try:
                result = minimize(
                    cost, 
                    l_initial, 
                    constraints=[monotonicity_constraint, target_constraint, boundary_constraint], 
                    tol=1e-8, 
                    options={'maxiter': 10000, 'disp': False}
                )
                
                if result.success:
                    target_costs.append(result.fun)
                else:
                    target_costs.append(np.nan)
            except:
                # Handle optimization failures gracefully
                target_costs.append(np.nan)
        else:
            # Target ratio is outside feasible bounds (only relevant when denominator is not reversible)
            target_costs.append(np.nan)
    else:
        target_costs.append(np.nan)

#=====================================
#7. Return results to Stata
#=====================================

#-------------------------------------
#7.1 Set local macros for each variable
#-------------------------------------

for i in range(num_vars):
    var_num = i + 1
    
    Macro.setLocal(f"orig_ratio_{var_num}", str(original_ratios[i]))
    Macro.setLocal(f"min_ratio_{var_num}", str(min_ratios[i]))
    Macro.setLocal(f"max_ratio_{var_num}", str(max_ratios[i]))
    
    if not np.isnan(target_costs[i]):
        Macro.setLocal(f"target_cost_{var_num}", str(target_costs[i]))
    else:
        Macro.setLocal(f"target_cost_{var_num}", ".")

#-------------------------------------
#7.2 Store summary information
#-------------------------------------

Macro.setLocal("num_variables", str(num_vars))

print(f"MRS reversal analysis completed for {num_vars} variables")