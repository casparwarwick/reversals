********************************************************************************
*Reversing the reversal
********************************************************************************
*mrs_reverser
********************************************************************************

cap program drop mrs_reverser
program mrs_reverser, rclass

	*=====================================
	*1. Initial setup
	*=====================================

	syntax,									/// 
	denom(varlist max=1)					/// Denominator variable (only one variable accepted)
	[										///
	python									/// Use Python numerical cost minimization instead of the analytical active-set search (default)
	target_ratio(real -999)					/// Target ratio for cost calculation. Only relevant when python is not specified, or in combination with alpha/theil.
	alpha(real 2)							/// Alpha parameter for cost function (default: 2). Only relevant when python is specified.
	theil									/// Use normalized Theil index as cost function (overrides alpha option). Only relevant when python is specified.
	keep(string) 							/// Specifies list of variables to be kept in the displayed results table
	MAX_attempts(integer 50000)				/// Cap on number of cra_solve_face calls in analytical search. Ignored when python is specified.
	]		

	qui {
	
	*-------------------------------------
	*	1.1 Option logic and warnings
	*-------------------------------------
	
	* Default is the analytical search. The `python' option switches to numerical
	* optimisation. Warn about options that are ignored under the chosen routine.
	
	if "`python'" == "" {
		* Analytical routine in use - python-only options are ignored.
		local _ignored_py ""
		if `alpha'   != 2  local _ignored_py "`_ignored_py' alpha()"
		if "`theil'" != "" local _ignored_py "`_ignored_py' theil"
		if "`_ignored_py'" != "" {
			noi dis as err "Note: the following option(s) are ignored when the python option is not specified:`_ignored_py'."
		}
	}
	else {
		* Python routine in use - analytical-only options are ignored.
		if `max_attempts' != 50000 {
			noi dis as err "Note: the max_attempts() option is ignored when the python option is specified."
		}
	}
	
	* If python is requested, check that numpy and scipy are available.
	* If not, fall back to the analytical routine.
	if "`python'" != "" {
		capture python which numpy
		if _rc != 0 {
			noi dis as err "Numpy not available. Falling back to the analytical search."
			local python ""
		}
	}
	if "`python'" != "" {
		capture python which scipy
		if _rc != 0 {
			noi dis as err "SciPy not available. Falling back to the analytical search."
			local python ""
		}
	}
	
	* Clean up any leftover matrices from previous runs
	cap mat drop _labels_depvar
	cap mat drop _numerator_coeffs
	cap mat drop _denominator_coeffs
	
	*-------------------------------------
	*	1.2 Housekeeping
	*-------------------------------------
		
	* Store the previous model
	tempname prevmodel
	estimates store `prevmodel'
	
	* Instantiate some local names to be used later
	tempname tmp_mat 		// temporary matrix used in various parts of the program
	tempname tmp_mat2 		// temporary matrix used in various parts of the program	
	tempname tmp_w_mat 		// records the coefficients from regressions of hd
	tempname d_result		// records signs of coefficients in regressions of hd
	
	*-------------------------------------
	*	1.3 Get original quantities needed later
	*-------------------------------------
	
	* Get the original estimates
	tempname orig_b
	matrix `orig_b' = e(b)
	tempname orig_b_backup
	matrix `orig_b_backup' = e(b)
			
	* Get the full original command and the dependent variable
	local full_command  "`e(cmdline)'"
	local depvar "`e(depvar)'"
	
	* Get a rescaled version of the original variable
	sum `depvar', meanonly
	local min_depvar = r(min)
	local max_depvar = r(max)
	
	local scale_min = `min_depvar'
	local scale_max = `max_depvar'
	local scale_range = `scale_max' - `scale_min'
			
	* Get the labels into a matrix (needed for the python routine)
	levelsof `depvar', matrow(_labels_depvar)
	
	*=====================================
	*2. Assess if reversals are possible for denominator
	*=====================================
	
	*-------------------------------------
	*2.1 Get hd 
	*-------------------------------------
	
	levelsof `depvar', local(levels_depvar)
	local n = 1 
	foreach i of local levels_depvar {
		tempvar dstub`n'
		gen 	`dstub`n'' = 0 if `depvar' != . & e(sample)==1
		qui replace `dstub`n'' = 1 if `depvar' <= `i'	& e(sample)==1
		local ++n
	}
	local n_levels_depvar = `n' -1

	*-------------------------------------
	*2.2 Run regressions of hd 
	*-------------------------------------

	local nrows_d_result = `n' - 2
	forvalues n=1(1)`nrows_d_result' {
		
		* Get the estimation command
		local to_run = subinword("`full_command'", "`depvar'", "`dstub`n''", 1)
		
		* Run the regression
		`to_run' 
		
		* Record the result
		matrix `tmp_mat' = e(b) 
		if `n'==1 matrix `tmp_w_mat'=`tmp_mat' 
		else	  matrix `tmp_w_mat' = (`tmp_w_mat' \ `tmp_mat')
		
		* Just get the signs and place in matrix
		mata : st_matrix("`tmp_mat2'", st_matrix("`tmp_mat'") :> 0)		
		if `n'==1 matrix `d_result'=`tmp_mat2'
		else	  matrix `d_result' = (`d_result' \ `tmp_mat2')
	
	}
		
	*-------------------------------------
	*2.3 Check whether signs of denom stay the same
	*-------------------------------------

	matrix colnames `d_result'=`:colnames `tmp_w_mat'' 
	tempname signs_denominator
	matselrc `d_result' `signs_denominator', c(`denom')
	
	tempname is_reversible_mat
	matrix `is_reversible_mat' = J(1,rowsof(`signs_denominator'),1)*`signs_denominator'
	
	* Check if denominator is reversible (set flag for later use)
	local denom_reversible = 0
	if `is_reversible_mat'[1,1] != 0 & `is_reversible_mat'[1,1] != `nrows_d_result' {
		local denom_reversible = 1
	}
		
	*=====================================
	*3. Prepare coefficient vectors for ratio analysis
	*=====================================

	*-------------------------------------
	*3.1 Get numerator variables (all except denominator and constant)
	*-------------------------------------
			
	local explanatory_vars "`:colnames `orig_b''"
	local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
	local numerator_vars = subinstr("`explanatory_vars'", "`denom'", "",1)
	
	* Clean up extra spaces
	local numerator_vars = trim("`numerator_vars'")
	local numerator_vars = subinstr("`numerator_vars'", "  ", " ",.)
	
	*-------------------------------------
	*3.2 Extract coefficients
	*-------------------------------------
	
	* Extract denominator coefficients (column vector of bd's for the denominator)
	matselrc `tmp_w_mat' _denominator_coeffs, c(`denom')
	
	* Extract numerator coefficients for each variable
	local var_counter = 0
	foreach var of local numerator_vars {
		local ++var_counter
		tempname num_`var_counter'
		matselrc `tmp_w_mat' `num_`var_counter'', c(`var')
		if `var_counter' == 1 {
			matrix _numerator_coeffs = `num_`var_counter''
		}
		else {
			matrix _numerator_coeffs = (_numerator_coeffs, `num_`var_counter'')
		}
	}
	
	*-------------------------------------
	*3.3 Set up locals
	*-------------------------------------
	
	* Target ratio (if specified)
	if `target_ratio' != -999 {
		local has_target_ratio = 1
		local target_ratio_value = `target_ratio'
	}
	else {
		local has_target_ratio = 0
		local target_ratio_value = 0
	}
	
	* Pass reversible flag to Python (used by Python script)
	local denom_reversible_flag = `denom_reversible'
	
	* Get original (linear) coefficients for the denominator and each numerator
	local orig_denom = `orig_b'[1, colnumb(`orig_b', "`denom'")]
	
	*=====================================
	*4. Run analysis
	*=====================================
	
	if "`python'" != "" {
		*-------------------------------------
		*4.1 Run Python optimization
		*-------------------------------------
		
		python script "`c(sysdir_plus)'py/mrs_reverser_python.py"
		
	}
	else {
		*-------------------------------------
		*4.2 Run analytical active-set search
		*-------------------------------------
		
		* For each numerator variable, compute:
		*   - orig_ratio_X      : original linear ratio (from orig_b)
		*   - min_ratio_X, max_ratio_X : bounds on the ratio across positive
		*                                 monotonic transformations
		*   - target_cost_X     : minimum SD-based cost to attain target_ratio_value,
		*                          conditional on it being feasible
		
		local n_gaps = `n_levels_depvar' - 1
		
		* Vector of dichotomised coefficients for the denominator
		tempname bd_denom_col
		matrix `bd_denom_col' = _denominator_coeffs
		
		local var_counter = 0
		foreach var of local numerator_vars {
			local ++var_counter
			
			* Original linear ratio
			local orig_num = `orig_b'[1, colnumb(`orig_b', "`var'")]
			if abs(`orig_denom') > 1e-15 {
				local orig_ratio_`var_counter' = `orig_num' / `orig_denom'
			}
			else {
				local orig_ratio_`var_counter' = .
			}
			
			* Pull the column of dichotomised numerator coefficients
			tempname bd_num_col
			matrix `bd_num_col' = _numerator_coeffs[1...,`var_counter']
			
			* Compute min and max ratios across k = 1..K-1
			mata: st_local("min_ratio_`var_counter'", strofreal(min(st_matrix("`bd_num_col'") :/ st_matrix("`bd_denom_col'"))))
			mata: st_local("max_ratio_`var_counter'", strofreal(max(st_matrix("`bd_num_col'") :/ st_matrix("`bd_denom_col'"))))
			
			* If the denominator is reversible (sign-flips across k), the implied
			* ratio is unbounded - report as +/- infinity.
			if `denom_reversible' == 1 {
				local min_ratio_`var_counter' = "-inf"
				local max_ratio_`var_counter' = "inf"
			}
			
			* Minimum cost to attain target_ratio_value
			if `has_target_ratio' == 1 {
				
				* The cost of attaining ratio rho is the cost of "sign-reversing"
				* the linear combination c_k = bd_num_k - rho * bd_denom_k to 0.
				* So we form c_k and feed it to cra_search.
				tempname c_vec
				mata: st_matrix("`c_vec'", st_matrix("`bd_num_col'") :- `target_ratio_value' :* st_matrix("`bd_denom_col'"))
				
				* Quick feasibility check: the target ratio must lie within the
				* (numeric) bounds on the achievable ratio.
				local feasible = 1
				if "`min_ratio_`var_counter''" != "-inf" {
					if `target_ratio_value' < `min_ratio_`var_counter'' local feasible = 0
				}
				if "`max_ratio_`var_counter''" != "inf" {
					if `target_ratio_value' > `max_ratio_`var_counter'' local feasible = 0
				}
				
				if `feasible' == 1 {
					tempname result_row
					mata: st_matrix("`result_row'", cra_search(st_matrix("`c_vec'"), `n_gaps', `max_attempts'))
					
					* Unpack the cost (first element of the result row)
					local target_cost_`var_counter' = `result_row'[1, 1]
					local n_attempts_`var_counter'  = `result_row'[1, `n_gaps' + 3]
					local hit_limit_`var_counter'   = `result_row'[1, `n_gaps' + 4]
				}
				else {
					local target_cost_`var_counter' = .
					local n_attempts_`var_counter'  = .
					local hit_limit_`var_counter'   = 0
				}
			}
			else {
				local target_cost_`var_counter' = .
				local n_attempts_`var_counter'  = .
				local hit_limit_`var_counter'   = 0
			}
		}
	}
	
	*=====================================
	*5. Display results
	*=====================================
	
	*-------------------------------------
	*5.1 Display header
	*-------------------------------------
	
	noi dis ""
	noi dis as text "{hline 78}"
	noi dis as text "Coefficient Ratios Relative to " as result "`denom'"
	noi dis as text "{hline 78}"
	noi dis ""
	
	*-------------------------------------
	*5.2 Determine which variables to display (respect keep() option)
	*-------------------------------------
	
	if "`keep'" != "" {
		local display_vars ""
		foreach v of local keep {
			* Only include vars that are actually in the numerator list
			local in_list : list v in numerator_vars
			if `in_list' local display_vars "`display_vars' `v'"
		}
		local display_vars = trim("`display_vars'")
		if "`display_vars'" == "" {
			noi dis as err "Note: none of the variables in keep() match numerator variables. Displaying all numerator variables."
			local display_vars "`numerator_vars'"
		}
	}
	else {
		local display_vars "`numerator_vars'"
	}
	
	*-------------------------------------
	*5.3 Display results table
	*-------------------------------------
	
	noi dis as text %15s "Variable" %15s "Orig.ratio" %15s "Min.ratio" %15s "Max.ratio" %15s "Min.cost"
	noi dis as text "{hline 78}"
	
	* Loop through numerator variables and display results. We always loop over the
	* full list of numerator vars internally to keep var_counter aligned with the
	* matrices, and only print rows for vars in display_vars.
	local var_counter = 0
	foreach var of local numerator_vars {
		local ++var_counter
		
		* Skip if not in display list
		local in_display : list var in display_vars
		if `in_display' == 0 continue
		
		* Get the ratios for this variable
		local orig_val = "`orig_ratio_`var_counter''"
		local min_val = "`min_ratio_`var_counter''"  
		local max_val = "`max_ratio_`var_counter''"
		local cost_val = "`target_cost_`var_counter''"
		
		* Format values
		if "`orig_val'" != "" & "`orig_val'" != "." {
			local orig_disp : di %9.3f `orig_val'
		}
		else local orig_disp = "."
		
		if "`min_val'" != "" {
			if "`min_val'" == "-inf" {
				local min_disp = "  -infty"
			}
			else {
				local min_disp : di %9.3f `min_val'
			}
		}  
		else local min_disp = "."
		
		if "`max_val'" != "" {
			if "`max_val'" == "inf" {
				local max_disp = "   infty"
			}
			else {
				local max_disp : di %9.3f `max_val'
			}
		}
		else local max_disp = "."
		
		if "`cost_val'" != "" & `cost_val' != . {
			local cost_disp : di %9.3f `cost_val'
		}
		else local cost_disp = "."
		
		* Display the row
		noi dis as result %15s "`var'" %15s "`orig_disp'" %15s "`min_disp'" %15s "`max_disp'" %15s "`cost_disp'"
	}
	
	noi dis as text "{hline 78}"
	noi dis ""
	
	*-------------------------------------
	*5.4 Display footer info
	*-------------------------------------
	
	if `has_target_ratio' == 1 {
		noi dis as text "Target ratio: " as result `target_ratio_value'
		noi dis ""
	}
	else {
		noi dis as text "Target ratio: " as result "no target ratio specified"
		noi dis ""
	}
	
	* Display warning if denominator is reversible
	if `denom_reversible' == 1 {
		noi dis as text "Warning: The specified denominator is reversible. Ratios are unbounded."
		noi dis ""
	}
	
	* If analytical search hit the cap anywhere, warn
	if "`python'" == "" & `has_target_ratio' == 1 {
		local any_hit = 0
		local var_counter = 0
		foreach var of local numerator_vars {
			local ++var_counter
			if "`hit_limit_`var_counter''" == "1" local any_hit = 1
		}
		if `any_hit' == 1 {
			noi dis as err "WARNING: For at least one variable, the analytical search reached the cap of `max_attempts' calls."
			noi dis as err "         Reported minimum costs may not be the global optimum for those variables."
			noi dis as err "         See r(hit_limit) for which variables are affected. You may want to increase max_attempts()."
			noi dis ""
		}
	}
	
	*=====================================
	*6. Store results in r()
	*=====================================
	
	*-------------------------------------
	*6.1 Create result matrices
	*-------------------------------------
	
	* Count numerator variables for matrix dimensions
	local var_count : word count `numerator_vars'
	
	* Create matrices for storing results
	tempname result_matrix ratio_matrix minratio_matrix maxratio_matrix
	tempname cost_matrix n_attempts_matrix hit_limit_matrix
	
	* Initialize matrices
	matrix `ratio_matrix' = J(`var_count', 1, .)
	matrix `minratio_matrix' = J(`var_count', 1, .)
	matrix `maxratio_matrix' = J(`var_count', 1, .)
	
	if `has_target_ratio' == 1 {
		matrix `cost_matrix' = J(`var_count', 1, .)
	}
	
	* Analytical-search-only diagnostics
	if "`python'" == "" {
		matrix `n_attempts_matrix' = J(`var_count', 1, .)
		matrix `hit_limit_matrix' = J(`var_count', 1, .)
	}
	
	* Create result display matrix
	if `has_target_ratio' == 1 {
		matrix `result_matrix' = J(`var_count', 4, .)
	}
	else {
		matrix `result_matrix' = J(`var_count', 3, .)
	}
	
	* Fill matrices with computed results
	local var_counter = 0
	foreach var of local numerator_vars {
		local ++var_counter
		
		* Store basic ratios
		if "`orig_ratio_`var_counter''" != "" & "`orig_ratio_`var_counter''" != "." {
			matrix `ratio_matrix'[`var_counter', 1] = `orig_ratio_`var_counter''
			matrix `result_matrix'[`var_counter', 1] = `orig_ratio_`var_counter''
		}
		
		* Handle infinite bounds
		if "`min_ratio_`var_counter''" == "-inf" {
			matrix `minratio_matrix'[`var_counter', 1] = -999999999
			matrix `result_matrix'[`var_counter', 2] = -999999999
		}
		else if "`min_ratio_`var_counter''" != "" {
			matrix `minratio_matrix'[`var_counter', 1] = `min_ratio_`var_counter''
			matrix `result_matrix'[`var_counter', 2] = `min_ratio_`var_counter''
		}
		
		if "`max_ratio_`var_counter''" == "inf" {
			matrix `maxratio_matrix'[`var_counter', 1] = 999999999
			matrix `result_matrix'[`var_counter', 3] = 999999999
		}
		else if "`max_ratio_`var_counter''" != "" {
			matrix `maxratio_matrix'[`var_counter', 1] = `max_ratio_`var_counter''
			matrix `result_matrix'[`var_counter', 3] = `max_ratio_`var_counter''
		}
		
		* Store target costs if applicable
		if `has_target_ratio' == 1 {
			if "`target_cost_`var_counter''" != "" & `target_cost_`var_counter'' != . {
				matrix `cost_matrix'[`var_counter', 1] = `target_cost_`var_counter''
				matrix `result_matrix'[`var_counter', 4] = `target_cost_`var_counter''
			}
		}
		
		* Analytical-search-only diagnostics
		if "`python'" == "" {
			if "`n_attempts_`var_counter''" != "" & "`n_attempts_`var_counter''" != "." {
				matrix `n_attempts_matrix'[`var_counter', 1] = `n_attempts_`var_counter''
			}
			if "`hit_limit_`var_counter''" != "" {
				matrix `hit_limit_matrix'[`var_counter', 1] = `hit_limit_`var_counter''
			}
		}
	}
	
	* Set matrix row and column names
	matrix rownames `ratio_matrix' = `numerator_vars'
	matrix rownames `minratio_matrix' = `numerator_vars'
	matrix rownames `maxratio_matrix' = `numerator_vars'
	matrix rownames `result_matrix' = `numerator_vars'
	
	matrix colnames `ratio_matrix' = "orig_ratio"
	matrix colnames `minratio_matrix' = "min_ratio"
	matrix colnames `maxratio_matrix' = "max_ratio"
	
	if `has_target_ratio' == 1 {
		matrix rownames `cost_matrix' = `numerator_vars'
		matrix colnames `cost_matrix' = "cost"
		matrix colnames `result_matrix' = "orig_ratio" "min_ratio" "max_ratio" "cost"
	}
	else {
		matrix colnames `result_matrix' = "orig_ratio" "min_ratio" "max_ratio"
	}
	
	if "`python'" == "" {
		matrix rownames `n_attempts_matrix' = `numerator_vars'
		matrix rownames `hit_limit_matrix' = `numerator_vars'
		matrix colnames `n_attempts_matrix' = "n_attempts"
		matrix colnames `hit_limit_matrix' = "hit_limit"
	}
	
	*-------------------------------------
	*6.2 Return results in r()
	*-------------------------------------
	
	* Main result matrices
	return matrix result = `result_matrix'
	return matrix ratio = `ratio_matrix'
	return matrix minratio = `minratio_matrix'
	return matrix maxratio = `maxratio_matrix'
	
	* Cost matrices (if target ratio specified)
	if `has_target_ratio' == 1 {
		return matrix cost = `cost_matrix'
	}
	
	* Analytical-search-only diagnostics
	if "`python'" == "" {
		return matrix n_attempts = `n_attempts_matrix'
		return matrix hit_limit = `hit_limit_matrix'
	}
	
	*=====================================
	*7. Clean up and restore
	*=====================================
	
	*-------------------------------------
	*7.1 Clean up matrices
	*-------------------------------------
	
	cap mat drop _labels_depvar
	cap mat drop _denominator_coeffs
	cap mat drop _numerator_coeffs

	*-------------------------------------
	*7.2 Restore the original model
	*-------------------------------------
	
	qui estimates restore `prevmodel'
	
	}	// ends the qui condition
	
end


*========================================
* Mata functions for the analytical search.
* These are identical to the routines used by coeff_reverser.
* Including them here too so that mrs_reverser can be used standalone.
*========================================

clear mata
mata:

real rowvector cra_solve_face(real colvector b, real scalar A_bits, real scalar n_gaps)
{
    // Closed-form solution of the reduced problem on a single face.
    // Returns (C, Delta_1, ..., Delta_n_gaps) as a row vector,
    // or a row of missing values if the reduced problem is degenerate
    // (M < 2 or zero variance of free dichotomised coefficients).

    real colvector active, free, Delta
    real scalar k, M, mu, V, C

    active = J(n_gaps, 1, 0)
    for (k = 1; k <= n_gaps; k++) {
        active[k] = mod(floor(A_bits / 2^(k - 1)), 2)
    }
    free = 1 :- active
    M    = sum(free)
    if (M < 2) return(J(1, n_gaps + 1, .))

    mu = sum(b :* free) / M
    V  = sum(((b :- mu):^2) :* free) / M
    if (V < 1e-15) return(J(1, n_gaps + 1, .))

    Delta = (1/M :- mu * (b :- mu) / (M * V)) :* free
    C     = sqrt(n_gaps / (n_gaps - 1) * sum((Delta :- 1/n_gaps):^2))
    return((C, Delta'))
}

real rowvector cra_naive_feasible(real colvector b, real scalar n_gaps)
{
    // Closed-form two-block warm start.

    real scalar S_pos, S_neg, n_pos, n_neg, n_zero, lambda, T, C, k
    real colvector Delta

    S_pos  = 0
    S_neg  = 0
    n_pos  = 0
    n_neg  = 0
    n_zero = 0
    for (k = 1; k <= n_gaps; k++) {
        if (b[k] > 0) {
            S_pos = S_pos + b[k]
            n_pos = n_pos + 1
        }
        else if (b[k] < 0) {
            S_neg = S_neg + b[k]
            n_neg = n_neg + 1
        }
        else {
            n_zero = n_zero + 1
        }
    }

    if (n_pos == 0 | n_neg == 0) return(J(1, n_gaps + 1, .))

    lambda = -S_pos / S_neg
    T      = lambda * n_neg + n_pos + n_zero

    Delta = J(n_gaps, 1, .)
    for (k = 1; k <= n_gaps; k++) {
        if (b[k] < 0) Delta[k] = lambda / T
        else          Delta[k] = 1 / T
    }

    C = abs(lambda - 1) / T * sqrt(n_neg * (n_pos + n_zero) / (n_gaps - 1))

    return((C, Delta'))
}

real rowvector cra_search(real colvector b, real scalar n_gaps, real scalar max_attempts)
{
    // Active-set search with cost-based pruning and a closed-form warm start.
    //
    // Returns 1 x (n_gaps + 4):
    //   (best_C, best_Delta_1, ..., best_Delta_n_gaps, best_A_bits,
    //    n_attempts, hit_limit).

    real scalar best_C, best_A_bits, max_level, level
    real scalar i, k, A, parent_A, child_A
    real scalar attempts, hit_limit
    real colvector best_Delta, Delta
    real colvector current_frontier, current_costs
    real colvector new_frontier, new_costs
    real colvector candidates, keep
    real rowvector sol, naive_sol

    max_level   = n_gaps - 1
    best_C      = .
    best_Delta  = J(n_gaps, 1, .)
    best_A_bits = -1
    attempts    = 0
    hit_limit   = 0

    // Level 0
    sol      = cra_solve_face(b, 0, n_gaps)
    attempts = attempts + 1
    if (sol[1, 1] == .) {
        return((., J(1, n_gaps, .), -1, attempts, hit_limit))
    }

    Delta = sol[1, 2..(n_gaps + 1)]'
    if (min(Delta) >= -1e-10) {
        return((sol[1, 1], Delta', 0, attempts, hit_limit))
    }

    // Closed-form naive warm start (free)
    naive_sol = cra_naive_feasible(b, n_gaps)
    if (naive_sol[1, 1] == .) {
        return((., J(1, n_gaps, .), -1, attempts, hit_limit))
    }
    best_C      = naive_sol[1, 1]
    best_Delta  = naive_sol[1, 2..(n_gaps + 1)]'
    best_A_bits = -2

    if (attempts >= max_attempts) {
        hit_limit = 1
        return((best_C, best_Delta', best_A_bits, attempts, hit_limit))
    }

    // BFS over faces
    current_frontier = J(1, 1, 0)
    current_costs    = J(1, 1, sol[1, 1])

    for (level = 1; level <= max_level; level++) {

        keep             = current_costs :< best_C
        current_frontier = select(current_frontier, keep)
        current_costs    = select(current_costs,    keep)
        if (rows(current_frontier) == 0) break

        candidates = J(0, 1, .)
        for (i = 1; i <= rows(current_frontier); i++) {
            parent_A = current_frontier[i]
            for (k = 1; k <= n_gaps; k++) {
                if (mod(floor(parent_A / 2^(k - 1)), 2) == 0) {
                    child_A    = parent_A + 2^(k - 1)
                    candidates = candidates \ child_A
                }
            }
        }
        if (rows(candidates) == 0) break
        candidates = uniqrows(candidates)

        new_frontier = J(0, 1, .)
        new_costs    = J(0, 1, .)
        for (i = 1; i <= rows(candidates); i++) {

            if (attempts >= max_attempts) {
                hit_limit = 1
                break
            }

            A        = candidates[i]
            sol      = cra_solve_face(b, A, n_gaps)
            attempts = attempts + 1
            if (sol[1, 1] == .) continue
            if (sol[1, 1] >= best_C) continue

            Delta = sol[1, 2..(n_gaps + 1)]'
            if (min(Delta) >= -1e-10) {
                best_C      = sol[1, 1]
                best_Delta  = Delta
                best_A_bits = A
            }
            else {
                new_frontier = new_frontier \ A
                new_costs    = new_costs    \ sol[1, 1]
            }
        }

        if (hit_limit) break

        current_frontier = new_frontier
        current_costs    = new_costs
    }

    return((best_C, best_Delta', best_A_bits, attempts, hit_limit))
}
end
