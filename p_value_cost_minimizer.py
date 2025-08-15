#*******************************************************************************
#Reversing the reversal
#*******************************************************************************
#Python routine to find p-value bounds using alternative approach
#Based on p_val_minimiser_v5.py but integrated with coeff_reverser cost function
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
from scipy import stats
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint, BFGS  

#=====================================
#2. Define cost function (from sign_reversal_cost_minimizer.py)
#=====================================

def cost(l):
    alpha_value = float(Macro.getLocal('alpha'))
    use_theil = Macro.getLocal('theil') != ''
    
    if use_theil:
        # Use normalized Theil index
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
        K = len(l)
        minl = l[0]
        maxl = l[K-1]
        
        dl = [0]*(K-1)
        for i in range(0,K-1):
            dl[i] = l[i+1] - l[i]
        maxdl = maxl - minl
        N = K - 1
        maxvar = (1/N - 1/N**2)*maxdl**2
        Edl = maxdl/N

        var_comp = [0]*(N)
        for i in range(0,N):
            var_comp[i] = (dl[i] - Edl)**2 
    
        var = (1/N)*np.sum(var_comp)
        
        # Apply alpha parameter: cost = (var/maxvar)^(1/alpha)
        exponent = 1/alpha_value
        cost = (var/maxvar)**exponent
            
        return cost

#=====================================
#3. Define functions to compute p-values
#=====================================

#-------------------------------------
#3.1. Define function to calculate beta from bd coefficients
#-------------------------------------

def calculate_betas(bds, labels):
    """Calculate coefficient from bd regression coefficients and labels"""
    for i in range(0,len(labels)-1):
        if i>0:
            j=i+1
            b_from_d = b_from_d + (labels[i]-labels[j])*bds[:,i][np.newaxis, :]
        else:
            b_from_d = (labels[0]-labels[1])*bds[:,0][np.newaxis, :]
    return b_from_d

#-------------------------------------
#3.2. Define function to calculate variance-covariance matrix
#-------------------------------------

def calculate_variance_covariance(n, k, X, eds, se_type, labels, W, W_vec):
    """Calculate variance-covariance matrix for different SE types"""
    
    # Calculate e_from_d (residuals)
    for i in range(0,len(labels)-1):
        if i>0:
            j = i+1
            e_from_d = e_from_d + (labels[i]-labels[j])*eds[:,i]
        else:
            e_from_d = (labels[0]-labels[1])*eds[:,0]

    # Basic standard errors (OLS)
    if se_type == 1:
        scalar = ((e_from_d @ W @ e_from_d)/(n-k))  # Since W_vec.sum() = n.
        varcov = scalar * np.linalg.inv(X.T @ W @ X)

    # Heteroskedasticity-robust standard errors
    elif se_type == 2:
        e_sq = e_from_d ** 2
        finite_sample_correction = n / (n - k)
        # Correct sandwich estimator: (X'WX)^(-1) * X'W * diag(e_sq) * W * X * (X'WX)^(-1)
        XtWX_inv = np.linalg.inv(X.T @ W @ X)
        # Create diagonal matrix from squared residuals weighted appropriately
        middle_term = X.T @ W @ np.diag(e_sq) @ W @ X
        varcov = finite_sample_correction * XtWX_inv @ middle_term @ XtWX_inv

    return np.diag(varcov)

#-------------------------------------
#3.3. Define function to calculate p-values
#-------------------------------------

def calculate_p_values(bds, labels_transformed, n, k, X, eds, se_type, df, W, W_vec):
    """Calculate p-values for given labels transformation"""
    beta = calculate_betas(bds, labels_transformed)
    varcov = calculate_variance_covariance(n, k, X, eds, se_type, labels_transformed, W, W_vec)
    SEs = np.sqrt(varcov)
    t = beta/SEs
    p = 2 * stats.t.sf(abs(t), df)
    return p

#=====================================
#4. Import data from Stata
#=====================================

#-------------------------------------
#4.1 Import basic parameters
#-------------------------------------

# Get regression sample size 
n = int(Macro.getLocal('N'))

# Get SE type
se_name = Macro.getLocal('se_name')

se_type = 1
if se_name == "robust": 
    se_type = 2

# Target p-value
target_p = float(Macro.getLocal('critval'))

#-------------------------------------
#4.2 Import X matrix and other data
#-------------------------------------

# Independent variables (X matrix)
variables = Macro.getLocal('variables')
X = np.asarray(Data.get(variables))

# Weight variable (always exists, normalized so sum = N)
W_vec = np.asarray(Data.get("_weightvar")).flatten()
W = np.diag(W_vec)

# Residuals from d regressions (use _hd_residual_* pattern)
eds_vars = []
n_d_regressions = int(Macro.getLocal('nrows_d_result'))
for i in range(1, n_d_regressions + 1):
    eds_vars.append(f"_hd_residual_{i}")
eds_names = " ".join(eds_vars)
eds = np.asarray(Data.get(eds_names))

# Coefficients from d regressions
bds = np.asarray(Matrix.get("_bds"))
bds = bds.T

# Get actual dimensions from the data
# Before transpose: bds is (n_d_regressions x k)  
# After transpose: bds is (k x n_d_regressions)
k = bds.shape[0]  # Number of coefficients in each d regression (from bds matrix after transpose)
n_d_regressions = bds.shape[1]  # Number of d regressions

# Degrees of freedom - use actual k from bds
df = n - k

# Number of labels
nlabs = len(eds[0,:])+1

# Original labels
l_original = np.array(range(1,nlabs+1,1))

# Scale bounds
scale_min = float(Macro.getLocal('scale_min'))
scale_max = float(Macro.getLocal('scale_max'))

#=====================================
#5. Set constraints
#=====================================

#-------------------------------------
#5.1 Set constraint that labels are weakly increasing
#-------------------------------------

# Define monotonicity constraints
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
#5.2 Set constraint that labels need to be between scale_min and scale_max
#-------------------------------------

tmp1 = [0]*nlabs
tmp1[0]=1
tmp2 = [0]*nlabs
tmp2[nlabs-1]=1

boundary_constraint = LinearConstraint([tmp1,tmp2], [scale_min,scale_max], [scale_min,scale_max])

#=====================================
#6. Define helper functions for optimization
#=====================================

#-------------------------------------
#6.1 Define functions to return the p-values for a specific coefficient
#-------------------------------------

# For the lower bound (minimize p-value)
def p_one_arg_min(labels_transformed, coeff_idx):   
    p = calculate_p_values(bds, labels_transformed, n, k, X, eds, se_type, df, W, W_vec)
    return p[0][coeff_idx]

# For the upper bound (maximize p-value)
def p_one_arg_max(labels_transformed, coeff_idx):   
    p = calculate_p_values(bds, labels_transformed, n, k, X, eds, se_type, df, W, W_vec)
    return -p[0][coeff_idx]

#-------------------------------------
#6.2 Define wrapper functions for cost minimization with p-value constraints
#-------------------------------------

def minimize_wrapper_min(target_p_val, coeff_idx):
    """Find minimum cost such that p-value <= target_p_val"""
    
    def p_constraint(labels_transformed):
        return p_one_arg_min(labels_transformed, coeff_idx)
    
    ratio_constraint_nonlinear = NonlinearConstraint(p_constraint, -np.inf, target_p_val, jac='2-point', hess=BFGS()) 
    result = minimize(cost, l_original, constraints=[monotonicity_constraint, ratio_constraint_nonlinear, boundary_constraint], tol=1e-8, options = {'maxiter': 10000, 'disp': False})
    return result

def minimize_wrapper_max(target_p_val, coeff_idx):
    """Find minimum cost such that p-value >= target_p_val"""
    
    def p_constraint(labels_transformed):
        return p_one_arg_min(labels_transformed, coeff_idx)
    
    ratio_constraint_nonlinear = NonlinearConstraint(p_constraint, target_p_val, np.inf, jac='2-point', hess=BFGS()) 
    result = minimize(cost, l_original, constraints=[monotonicity_constraint, ratio_constraint_nonlinear, boundary_constraint], tol=1e-8, options = {'maxiter': 10000, 'disp': False})
    return result

#=====================================
#7. Import p-value bounds from Stata (already computed in section 2.5)
#=====================================

# Import the already computed p-value bounds from Stata
min_pval_matrix = np.asarray(Matrix.get("_min_pval"))
max_pval_matrix = np.asarray(Matrix.get("_max_pval"))

# Flatten to 1D arrays for consistency with rest of code
lower_final = min_pval_matrix.flatten()
upper_final = max_pval_matrix.flatten()

# Get original p-values for cost calculation
test_p = calculate_p_values(bds, l_original, n, k, X, eds, se_type, df, W, W_vec)
actual_k = len(test_p[0])  # Use actual number of p-values returned

#=====================================
#8. Find cost associated with reaching pre-specified p-value
#=====================================

# Find original p-value (reuse the test calculation)
orig_p = test_p

# Compute costs for each coefficient - use actual number of coefficients
costs = [0] * actual_k
for h in range(0, actual_k):
    # Check if p-value is within bounds
    if lower_final[h] <= target_p <= upper_final[h]:
        # Determine which optimizer to use based on original p-value
        if orig_p[0][h] > target_p:
            # Need to decrease p-value
            result = minimize_wrapper_min(target_p, h)
        else:
            # Need to increase p-value
            result = minimize_wrapper_max(target_p, h)
            
        costs[h] = result.fun
    else:
        # If target p-value is outside bounds, set to missing
        costs[h] = np.nan

#=====================================
#9. Store results back to Stata
#=====================================

Matrix.store("_orig", orig_p[0])
Matrix.store("_costs", costs)