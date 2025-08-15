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

	syntax, 								/// 
	denom(varlist max=1)					/// Denominator variable (only one variable accepted)
	[										///
	pythonno								/// Use routine with exponential transformations instead of Python cost minimization
	start(real -2) end(real 2) 				/// Specifies the smallest and largest value of c over which should be searched. Only relevant when pythonno is specified. 
	PRECision(real 0.1) 					/// Specifies the 'density' or 'precision' of the grid of c. Only relevant when pythonno is specified. 
	target_ratio(real -999)					/// Target ratio for cost calculation
	alpha(real 2)							/// Alpha parameter for cost function (default: 2)
	theil									/// Use normalized Theil index as cost function (overrides alpha option)
	keep(string) 							/// Specifies list of variables to be kept in the displayed results table
	]		

	qui {
	
	*-------------------------------------
	*	1.1 Set default behavior and option logic
	*-------------------------------------
	
	* Default: use Python for cost minimization, but check if it works. If not, revert to pythonno. 
	if "`pythonno'" == "" {
		capture python which numpy
		if _rc != 0 {
			noi dis "Numpy not available. Reverting to pythonno option."
			local pythonno "pythonno"
		}
		capture python which scipy
		if _rc != 0 {
			noi dis "SciPy not available. Reverting to pythonno option."
			local pythonno "pythonno"
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
	*4. Prepare matrices for Python analysis
	*=====================================

	*-------------------------------------
	*4.1 Get numerator variables (all except denominator and constant)
	*-------------------------------------
			
	local explanatory_vars "`:colnames `orig_b''"
	local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
	local numerator_vars = subinstr("`explanatory_vars'", "`denom'", "",1)
	
	* Clean up extra spaces
	local numerator_vars = trim("`numerator_vars'")
	local numerator_vars = subinstr("`numerator_vars'", "  ", " ",.)
	
	*-------------------------------------
	*4.2 Extract coefficients for Python
	*-------------------------------------
	
	* Extract denominator coefficients
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
	*4.3 Set up locals for Python
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
	
	* Pass reversible flag to Python
	local denom_reversible_flag = `denom_reversible'
	
	*=====================================
	*5. Run analysis (Python or exponential search)
	*=====================================
	
	if "`pythonno'" == "" {
		*-------------------------------------
		*5.1 Run Python optimization
		*-------------------------------------
		
		python script "`c(sysdir_plus)'py/mrs_reverser_python.py"
		
	}
	else {
		*-------------------------------------
		*5.2 Run exponential search routine
		*-------------------------------------
		
		* Store original coefficient bounds for later comparison
		local explanatory_vars "`:colnames `orig_b''"
		local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
		local numerator_vars = subinstr("`explanatory_vars'", "`denom'", "",1)
		local numerator_vars = trim("`numerator_vars'")
		local numerator_vars = subinstr("`numerator_vars'", "  ", " ",.)
		
		* Count numerator variables
		local num_vars : word count `numerator_vars'
		
		* Initialize result matrices
		tempname b_result_all
		tempname c_values
		
		*Tell the user how many regressions need to be run
		local N_c = ((`end' - `start')/`precision')+1
		noisily dis ""
		noisily dis "Total number of c-values to check: `N_c'. Progress:"

		*A counter to keep track of things
		local iter = 1
		
		*-------------------------------------
		*5.2.1 Loop through c values
		*-------------------------------------
		
		forvalues c=`start'(`precision')`end'{
			
			*-------------------------------------
			*5.2.2 Get minus sign before the exp() function for negative c
			*-------------------------------------

			if `c' < 0 	local minus "-"
			if `c' > 0	local minus ""
			
			*-------------------------------------
			*5.2.3 Create exponential transformation
			*-------------------------------------

			tempvar trans_depvar_rescaled
			if (`c' < 0.0000001 & `c' > -0.0000001) gen `trans_depvar_rescaled' = `depvar_rescaled' // I.e. if c==0
			else {
				gen `trans_depvar_rescaled' = ((`minus'exp(`depvar'*`c') - `minus'exp(`min_depvar'*`c')) / (`minus'exp(`max_depvar'*`c') - `minus'exp(`min_depvar'*`c')))*(`scale_max'-`scale_min') + `scale_min'
			}

			*-------------------------------------
			*5.2.4 Get a matrix to record differences in labels
			*-------------------------------------

			tempname dlabels
			matrix `dlabels' = J(`nrows_d_result',`nrows_d_result',0)
			local n = 0
			levelsof `trans_depvar_rescaled', local(levels_transdepvar_rescales)
			foreach l of local levels_transdepvar_rescales {
				if `n'==0 {
					local k = `l' 	// k here denotes the previous label. Do nothing in the first iteration and just update k for the next iteration.
				}
				else { 				// produces a diagonal matrix with the right dimensions to multiply.
					cap mat `dlabels'[`n',`n'] = `k'-`l'
					local k = `l' 	// Store l in k for the next iteration. 
				}
				local ++n
			}
			
			*-------------------------------------
			*5.2.5 Use the results from regressions of hd to arrive at the transformed coefficients
			*-------------------------------------

			*Take the dmat result and multiply each element by the matrix recording differences of labels
			tempname tmp_coeff_result
			matrix `tmp_coeff_result' = `dlabels' * `tmp_w_mat'
			
			*Take the column-wise sum:
			matrix `tmp_coeff_result' = J(1,rowsof(`tmp_coeff_result'),1)*`tmp_coeff_result' 
			
			*Store result in the b_result matrix
			if `c'==`start' {
				matrix `b_result_all' = (`tmp_coeff_result')
				matrix `c_values' = (`c')
			}
			else {
				matrix `b_result_all' = (`b_result_all' \ `tmp_coeff_result')
				matrix `c_values' = (`c_values' \ `c')
			}
			
			*Display progress
			noisily _dots `iter' 0
			local ++iter
			drop `trans_depvar_rescaled'
		}
		
		*-------------------------------------
		*5.2.6 Calculate ratios and find bounds for each numerator variable
		*-------------------------------------
		
		* Extract denominator coefficients from b_result_all
		matselrc `b_result_all' exp_denom_coeffs, c(`denom')
		
		* Calculate ratios for each numerator variable
		local var_counter = 0
		foreach var of local numerator_vars {
			local ++var_counter
			
			* Extract numerator coefficients for this variable
			matselrc `b_result_all' exp_num_coeffs, c(`var')
			
			* Calculate ratios
			tempname ratio_matrix
			mata : st_matrix("`ratio_matrix'", st_matrix("exp_num_coeffs") :/ st_matrix("exp_denom_coeffs"))
			
			* Find min and max ratios
			mata : st_numscalar("min_ratio_sc", min(st_matrix("`ratio_matrix'")))
			mata : st_numscalar("max_ratio_sc", max(st_matrix("`ratio_matrix'")))
			
			* Get original ratio (find c closest to 0 for linear case)
			local closest_to_zero_idx = 1
			local min_abs_c = abs(`c_values'[1,1])
			forvalues i = 2/`=rowsof(`c_values')' {
				local current_abs_c = abs(`c_values'[`i',1])
				if `current_abs_c' < `min_abs_c' {
					local min_abs_c = `current_abs_c'
					local closest_to_zero_idx = `i'
				}
			}
			local orig_ratio_val = `ratio_matrix'[`closest_to_zero_idx',1]
			
			* Store results in locals for display
			local orig_ratio_`var_counter' = `orig_ratio_val'
			local min_ratio_`var_counter' = scalar(min_ratio_sc)
			local max_ratio_`var_counter' = scalar(max_ratio_sc)
			
			* Calculate target cost if specified
			if `has_target_ratio' == 1 {
				* Find c value that gives target ratio (if any)
				tempname target_found
				local target_found = 0
				local target_c = .
				
				forvalues i = 1/`=rowsof(`ratio_matrix')' {
					local current_ratio = `ratio_matrix'[`i',1]
					local current_c = `c_values'[`i',1]
					
					* Check if this ratio is close to target
					if abs(`current_ratio' - `target_ratio_value') < 0.001 {
						local target_found = 1
						local target_c = abs(`current_c')
						continue, break
					}
				}
				
				if `target_found' == 1 {
					local target_cost_`var_counter' = `target_c'
				}
				else {
					local target_cost_`var_counter' = .
				}
			}
			else {
				local target_cost_`var_counter' = .
			}
		}
		
		* Set number of variables for display
		local num_variables = `num_vars'
	}
	
	*-------------------------------------
	*5.3 Get results (from Python or exponential search)
	*-------------------------------------
	
	if "`pythonno'" == "" {
		* Get basic ratios from Python
		local orig_ratio_1 "`orig_ratio_1'"
		local min_ratio_1 "`min_ratio_1'"
		local max_ratio_1 "`max_ratio_1'"
		
		* Additional variables if they exist
		if "`orig_ratio_2'" != "" {
			local orig_ratio_2 "`orig_ratio_2'"
			local min_ratio_2 "`min_ratio_2'"  
			local max_ratio_2 "`max_ratio_2'"
		}
		
		local num_variables "`num_variables'"
	}
	* For exponential search, results are already stored in locals from section 5.2.6
	
	*=====================================
	*6. Display results
	*=====================================
	
	*-------------------------------------
	*6.1 Display header
	*-------------------------------------
	
	noi dis ""
	noi dis as text "{hline 78}"
	noi dis as text "Coefficient Ratios Relative to " as result "`denom'"
	noi dis as text "{hline 78}"
	noi dis ""
	
	*-------------------------------------
	*6.2 Display results table
	*-------------------------------------
	
	if "`pythonno'"=="" {
		noi dis as text %12s "Variable" %15s "Orig.ratio" %15s "Min.ratio" %15s "Max.ratio" %15s "Min.cost"
	}
	else {
		noi dis as text %12s "Variable" %15s "Orig.ratio" %15s "Min.ratio" %15s "Max.ratio" %15s "Min.c"
	}
	noi dis as text "{hline 78}"
	
	* Loop through numerator variables and display results
	local var_counter = 0
	foreach var of local numerator_vars {
		local ++var_counter
		
		* Get the ratios for this variable
		local orig_val = "`orig_ratio_`var_counter''"
		local min_val = "`min_ratio_`var_counter''"  
		local max_val = "`max_ratio_`var_counter''"
		local cost_val = "`target_cost_`var_counter''"
		
		* Format values
		if "`orig_val'" != "" {
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
		noi dis as result %12s "`var'" %15s "`orig_disp'" %15s "`min_disp'" %15s "`max_disp'" %15s "`cost_disp'"
	}
	
	noi dis as text "{hline 78}"
	noi dis ""
	
	if `has_target_ratio' == 1 {
		noi dis as text "Target ratio: " as result `target_ratio_value'
		if "`pythonno'" == "" {
			noi dis as text "Cost function: " as result "`cost_function_type'"
			if "`cost_function_type'" == "alpha" {
				noi dis as text "Alpha parameter: " as result `alpha'
			}
		}
		else {
			noi dis as text "Transformation: " as result "exponential search f(y)=exp(c*y)"
			noi dis as text "Search range: " as result "c ∈ [`start', `end'], precision = `precision'"
			noi dis as text "Note: " as result "Missing values may indicate target ratio falls outside search range"
		}
		noi dis ""
	}
	else {
		noi dis as text "Target ratio: " as result "no target ratio specified"
		if "`pythonno'" == "" {
			noi dis as text "Cost function: " as result "not applicable"
		}
		else {
			noi dis as text "Transformation: " as result "exponential search f(y)=exp(c*y)"
			noi dis as text "Search range: " as result "c ∈ [`start', `end'], precision = `precision'"
		}
		noi dis ""		
	}
	
	* Display warning if denominator is reversible (at the end)
	if `denom_reversible' == 1 {
		noi dis as text "Warning: The specified denominator is reversible. Results may not be accurate."
		noi dis ""
	}
	
	*=====================================
	*7. Store results in r()
	*=====================================
	
	*-------------------------------------
	*7.1 Create result matrices
	*-------------------------------------
	
	* Count numerator variables for matrix dimensions
	local var_count : word count `numerator_vars'
	
	* Create matrices for storing results
	tempname result_matrix ratio_matrix minratio_matrix maxratio_matrix
	tempname cost_matrix minc_matrix
	
	* Initialize matrices
	matrix `ratio_matrix' = J(`var_count', 1, .)
	matrix `minratio_matrix' = J(`var_count', 1, .)
	matrix `maxratio_matrix' = J(`var_count', 1, .)
	
	if `has_target_ratio' == 1 {
		if "`pythonno'" == "" {
			matrix `cost_matrix' = J(`var_count', 1, .)
		}
		else {
			matrix `minc_matrix' = J(`var_count', 1, .)
		}
	}
	
	* Create result display matrix
	if `has_target_ratio' == 1 {
		matrix `result_matrix' = J(`var_count', 5, .)
	}
	else {
		matrix `result_matrix' = J(`var_count', 4, .)
	}
	
	* Fill matrices with computed results
	local var_counter = 0
	foreach var of local numerator_vars {
		local ++var_counter
		
		* Store basic ratios
		matrix `ratio_matrix'[`var_counter', 1] = `orig_ratio_`var_counter''
		
		* Handle infinite bounds
		if "`min_ratio_`var_counter''" == "-inf" {
			matrix `minratio_matrix'[`var_counter', 1] = -999999999
		}
		else {
			matrix `minratio_matrix'[`var_counter', 1] = `min_ratio_`var_counter''
		}
		
		if "`max_ratio_`var_counter''" == "inf" {
			matrix `maxratio_matrix'[`var_counter', 1] = 999999999
		}
		else {
			matrix `maxratio_matrix'[`var_counter', 1] = `max_ratio_`var_counter''
		}
		
		* Fill result display matrix
		matrix `result_matrix'[`var_counter', 1] = `orig_ratio_`var_counter''
		
		if "`min_ratio_`var_counter''" == "-inf" {
			matrix `result_matrix'[`var_counter', 2] = -999999999
		}
		else {
			matrix `result_matrix'[`var_counter', 2] = `min_ratio_`var_counter''
		}
		
		if "`max_ratio_`var_counter''" == "inf" {
			matrix `result_matrix'[`var_counter', 3] = 999999999
		}
		else {
			matrix `result_matrix'[`var_counter', 3] = `max_ratio_`var_counter''
		}
		
		* Store target costs if applicable
		if `has_target_ratio' == 1 {
			if "`pythonno'" == "" {
				if `target_cost_`var_counter'' != . {
					matrix `cost_matrix'[`var_counter', 1] = `target_cost_`var_counter''
					matrix `result_matrix'[`var_counter', 4] = `target_cost_`var_counter''
				}
			}
			else {
				if `target_cost_`var_counter'' != . {
					matrix `minc_matrix'[`var_counter', 1] = `target_cost_`var_counter''
					matrix `result_matrix'[`var_counter', 4] = `target_cost_`var_counter''
				}
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
		if "`pythonno'" == "" {
			matrix rownames `cost_matrix' = `numerator_vars'
			matrix colnames `cost_matrix' = "cost"
			matrix colnames `result_matrix' = "orig_ratio" "min_ratio" "max_ratio" "cost"
		}
		else {
			matrix rownames `minc_matrix' = `numerator_vars'
			matrix colnames `minc_matrix' = "min_c"
			matrix colnames `result_matrix' = "orig_ratio" "min_ratio" "max_ratio" "min_c"
		}
	}
	else {
		matrix colnames `result_matrix' = "orig_ratio" "min_ratio" "max_ratio"
	}
	
	*-------------------------------------
	*7.2 Return results in r()
	*-------------------------------------
	
	* Main result matrices
	return matrix result = `result_matrix'
	return matrix ratio = `ratio_matrix'
	return matrix minratio = `minratio_matrix'
	return matrix maxratio = `maxratio_matrix'
	
	* Cost matrices (if target ratio specified)
	if `has_target_ratio' == 1 {
		if "`pythonno'" == "" {
			return matrix cost = `cost_matrix'
		}
		else {
			return matrix minc = `minc_matrix'
		}
	}
	
	*=====================================
	*8. Clean up and restore
	*=====================================
	
	*-------------------------------------
	*8.1 Clean up matrices
	*-------------------------------------
	
	cap mat drop _labels_depvar
	cap mat drop _denominator_coeffs
	cap mat drop _numerator_coeffs

	*-------------------------------------
	*8.2 Restore the original model
	*-------------------------------------
	
	qui estimates restore `prevmodel'
	
	}	// ends the qui condition
	
end
