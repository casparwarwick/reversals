********************************************************************************
*Reversing the reversal
********************************************************************************
*coeff_reverser
********************************************************************************

cap program drop coeff_reverser
program coeff_reverser, rclass

	*=====================================
	*1. Initial setup
	*=====================================

	syntax, [ 								/// 
	pythonno								/// Use routine with exponential transformations instead of the Python cost minimization
	start(real -2) end(real 2) 				/// Specifies the smallest and largest value of c over which should be searched. Only relevant when pythonno is specified. 
	PRECision(real 0.1) 					/// Specifies the 'density' or 'precision' of the grid of c. For example precision(0.1) says that we evaluate values of c in steps of 0.1. Only relevant when pythonno is specified. 
	pvalue									/// display p-value statistics (min/max p-values, and minimum costs. If pythonno is not specified this only works with vanilla 'OLS regression')
	critval(real 0.05) 						/// Specifies the alpha level where we speak of statistical significance. Only relevant when pvalue is specified. 
	alpha(real 2)							/// Specifies the alpha parameter for the cost function (default: 2)
	theil									/// Use normalized Theil index as cost function (overrides alpha option)  
	revpoint(real 0)					/// Specifies the target value for sign reversal (default: 0)
	keep(string) 							/// Specifies list of variables to be kept in the displayed results table(s).
	dstub(string) 							/// Specifies that the binary dummy should be saved and storted in a stub specified by string.
	]		

	qui {
	
	*-------------------------------------
	*	1.1 Set default behavior and option logic
	*-------------------------------------
	
	* Default: use Python for cost minimization, but check if it works.
	if "`pythonno'" == "" {
		dis "test"
	}
	
	* IF NOT: Default: use fast routine (unless pythonno + pvalue specified)
	if "`pythonno'" != "" & "`pvalue'" != "" local fast ""
	else local fast "fast"
		
	*-------------------------------------
	*	1.2 Housekeeping
	*-------------------------------------
		
	*Store the previous model
	tempname prevmodel
	estimates store `prevmodel'
	
	*Instantiate some local names to be used later. 
	tempname b_result 		// records (transformed) coefficients
	tempname p_result 		// records p-values from transformed regressions
	tempname hd_p_vals		// records p-values from regressions of hd
	tempname tmp1 			// temporary matrix used in various parts of the program
	tempname tmp2			// temporary matrix used in various parts of the program
	tempname tmp_mat 		// temporary matrix used in various parts of the program
	tempname tmp_mat2		// temporary matrix used in various parts of the program
	tempname d_result 		// records which regressions of hd are responsible for reversal
	tempname dlabels 		// record differences in labels. Needed for "fast" routine
	tempname tmp_w_mat 		// records the coefficients from regressions of hd. 
	
	*-------------------------------------
	*	1.3 Get original quantities needed later
	*-------------------------------------
	
	*Get the original estimates
	tempname orig_b
	matrix `orig_b' = e(b)
	tempname orig_b_backup
	matrix `orig_b_backup' = e(b)
	
	tempname orig_v
	matrix `orig_v' = e(V)
	
	*Get N and number of parameters P
	local P = colsof(e(b)) 
	local N = e(N)
	
	*If the cluster option is specified record their number
	local ncluster = e(N_clust)
	
	*Get the full original command and the dependent variable
	local full_command  "`e(cmdline)'"
	local depvar "`e(depvar)'"
	
	*Get the vce type
	local se_name "`e(vcetype)'"
	if "`se_name'" == "Robust" local se_name = "robust"
	
	*Get the scale of the original variable
	sum `depvar', meanonly
	local scale_min = r(min)
	local scale_max = r(max)
					
	*Get the labels into a matrix (needed for the python routine)
	levelsof `depvar', matrow(_labels_depvar)
	
	*Get the signs of the original estimates
	tempname signs_mat signs_mat_pos signs_mat_neg
	mata : st_matrix("`signs_mat_pos'", st_matrix("`orig_b'") :> 0) 		// everything that is positive is 1
	mata: st_matrix("`signs_mat_neg'", (st_matrix("`orig_b'") :< 0) :* -1)  // everything that is negative is -1
	matrix `signs_mat' = `signs_mat_pos' + `signs_mat_neg' 					// adds the two to get positive and negative signs. Working with Stata matrices is tedious. 
		
	*=====================================
	*2. Assess if reversals are at all possible
	*=====================================
	
	*-------------------------------------
	*2.1 Get hd 
	*-------------------------------------
	
	levelsof `depvar', local(levels_depvar) // gets levels of depvar. Also used in later sections. 
	local n = 1 							// initialises counter (since depvar may not start at 1.)	
	foreach i of local levels_depvar {
		tempvar dstub`n'
		gen 	`dstub`n'' = 0 if `depvar' != . & e(sample)==1
		qui replace `dstub`n'' = 1 if `depvar' <= `i'	& e(sample)==1
		if "`dstub'" != "" gen `dstub'`n' =`dstub`n''
		local ++n
	}
	local n_levels_depvar = `n' -1

	*-------------------------------------
	*2.2 Run regressions of hd and check whether signs stay the same
	*-------------------------------------

	local nrows_d_result = `n' - 2 	// number of regressions to be run. -2 because the local n, defined above, increments one more time than is intuitive. 
	forvalues n=1(1)`nrows_d_result' {
		
		*Get the estimation command
		local to_run = subinword("`full_command'", "`depvar'", "`dstub`n''", 1) 	// changes the dependent variable of the command originally run.
		
		*Run the regression. 
		`to_run' 
		
		*Record the result
		matrix `tmp_mat' = e(b) 
		if `n'==1 matrix `tmp_w_mat'=`tmp_mat' 
		else	  matrix `tmp_w_mat' = (`tmp_w_mat' \ `tmp_mat')		

		*Just get the signs	and place in matrix
		mata : st_matrix("`tmp_mat2'", st_matrix("`tmp_mat'") :> 0)		
		if `n'==1 matrix `d_result'=`tmp_mat2'
		else	  matrix `d_result' = (`d_result' \ `tmp_mat2')
		
		*Save p-values
		tempname hd_p_vals_tmp
		tempname hd_covariance
		matrix `hd_covariance' = vecdiag(e(V))
		
		if `ncluster'!=. mata : st_matrix("`hd_p_vals_tmp'", 2*(J(1,`P',1)-t((`ncluster'-1), abs(st_matrix("`tmp_mat'") :/ sqrt(st_matrix("`hd_covariance'"))))))  		 
		if `ncluster'==. mata : st_matrix("`hd_p_vals_tmp'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`tmp_mat'") :/ sqrt(st_matrix("`hd_covariance'")))))) 
		
		if `n'==1 matrix `hd_p_vals'=`hd_p_vals_tmp' 
		else	  matrix `hd_p_vals' = (`hd_p_vals' \ `hd_p_vals_tmp')
		
		*Store residuals for p-value analysis (if needed)
		if "`pvalue'" != "" {
			cap drop _hd_residual_`n'
			predict _hd_residual_`n', residual
		}
	}
	
	*-------------------------------------
	*2.3 Get a matrix that simply indicates reversibility to revpoint (0=No; 1=Yes)
	*-------------------------------------
	
	*Check if revpoint is within bounds for each coefficient
	tempname is_reversible
	preserve
	clear
	svmat `tmp_w_mat', names("bd")
	ds
	local bd_vars = "`r(varlist)'"
	local n_vars = `: word count `bd_vars''
	
	*For each variable, check if revpoint is within [min, max] bounds
	foreach var in `bd_vars' {
		sum `var' if `var'!=., meanonly
		local min_val = r(min)
		local max_val = r(max)
		gen reversible_`var' = (`revpoint' <= -1*(`scale_max'-`scale_min')*`min_val' & `revpoint' >= -1*(`scale_max'-`scale_min')*`max_val') in 1
	}
	keep if _n == 1
	mkmat reversible_*, matrix(`is_reversible')
	restore 			

	*-------------------------------------
	*2.4 Get a separate matrix for zero-crossing reversibility (needed for p-value bounds)
	*-------------------------------------
	
	*Check if zero is within bounds for each coefficient
	tempname is_reversible_zero
	preserve
	clear
	svmat `tmp_w_mat', names("bd")
	ds
	local bd_vars = "`r(varlist)'"
	
	*For each variable, check if 0 is within [min, max] bounds
	foreach var in `bd_vars' {
		sum `var' if `var'!=., meanonly
		local min_val = r(min)
		local max_val = r(max)
		gen reversible_zero_`var' = (0 <= -1*(`scale_max'-`scale_min')*`min_val' & 0 >= -1*(`scale_max'-`scale_min')*`max_val') in 1
	}
	keep if _n == 1
	mkmat reversible_zero_*, matrix(`is_reversible_zero')
	restore

	*-------------------------------------
	*2.5 Get minimum and maximum p-value (only if pvalue option specified)
	*-------------------------------------
	
	if "`pvalue'" != "" {
		**Get maximum p-value from hd regressions
		preserve
		clear
		svmat `hd_p_vals', names("tmp")
		ds
		local vars = "`r(varlist)'"
		local n = 1
		foreach var in `vars' {
			sum `var'
			replace `var' = r(max)
			replace `var' = 1 if `is_reversible_zero'[1,`n'] == 1
			local ++n
		}
		keep if _n ==1
		tempname max_pval
		mkmat `vars', matrix(`max_pval')
		restore

		**Get minimum p-value from hd regressions
		preserve
		clear
		svmat `hd_p_vals', names("tmp")
		ds
		local vars = "`r(varlist)'"
		foreach var in `vars' {
			sum `var'
			replace `var' = r(min) 	
		}
		keep if _n ==1
		tempname min_pval
		mkmat `vars', matrix(`min_pval')
		restore
	}	
	
	*=====================================
	*3. Implement the cost-function approach 
	*=====================================

	if "`pythonno'"=="" {

		*-------------------------------------
		*3.1 Some setups
		*-------------------------------------
				
		local explanatory_vars "`:colnames `orig_b''"
		local explanatory_vars = subinstr("`explanatory_vars'", "_cons", "",1)
		python clear
		
		*-------------------------------------
		*3.2 Loop over explnatory variable
		*-------------------------------------
		
		preserve
		clear 
		svmat `tmp_w_mat', names("bd")
 
		set obs `n_levels_depvar'
			
		local n = 1
		foreach explanatory_var of local explanatory_vars {
			
			*-------------------------------------
			*3.3 Prepare things for python
			*-------------------------------------

			*Skip variable if no reversal can be achieved (i.e. when revpoint is outside the bounds of bd coefficients)
			sum bd`n' if bd`n'!=., meanonly
			local min_bd = r(min)
			local max_bd = r(max)
			
			*Check if revpoint is within the bounds [min_bd, max_bd]
			if (`revpoint' > -1*(`scale_max'-`scale_min')*`min_bd') | (`revpoint' < -1*(`scale_max'-`scale_min')*`max_bd') {
				gen new_labels`n' = .
				gen cost`n' = .
				local ++n
				continue
			}
			
			*Make sure that these variables, which Python will create, don't already exist. 
			cap drop python_labels
			cap drop python_cost
			
			*Get a sign variable to be imported into Python.
			cap drop sign
			gen sign = `signs_mat'[1,`n']
			
			*Get a variable with the coefficients on hd to be imported into Python.
			cap drop bd
			gen bd = bd`n'
			
			*Get reversal point to be imported into Python.
			cap drop revpoint
			gen revpoint = `revpoint'

			*-------------------------------------
			*3.4 Do the substantive Python bits
			*-------------------------------------
			
			noi python script "`c(sysdir_plus)'py/sign_reversal_cost_minimizer.py"
			
			*-------------------------------------
			*3.5 Get results into the right variables
			*-------------------------------------
			
			gen new_labels`n' = python_labels
			gen cost`n' = python_cost
			
			*-------------------------------------
			*3.6 Iterate counter
			*-------------------------------------

			local ++n		
		}
		
		*-------------------------------------
		*3.7 Save results to matrix
		*-------------------------------------
		
		local m = `n' -1
		local costs
		forvalues j = 1(1)`m' {
			local costs `costs' cost`j'
		}
		tempname cost_mat
		mkmat `costs', matrix(`cost_mat')
		matrix `cost_mat' = `cost_mat'[1,1...]
		
		local labels
		forvalues j = 1(1)`m' {
			local labels `labels' new_labels`j'
		}
		tempname label_mat
		mkmat `labels', matrix(`label_mat')
		
		restore 
		
		*=====================================
		*3.8 Compute p-value bounds (if pvalue option specified)
		*=====================================
		
		if "`pvalue'" != "" {
			
			*Store the current dataset state
			preserve
			
			*Keep only the regression sample
			keep if e(sample)==1
			
			*Get the weighting variable (if it exists, otherwise=1)
			cap drop _weightvar
			gen _weightvar = 1
			local weightvar = e(wexp)
			local weightvar = substr("`weightvar'", 3,.)
			if "`weightvar'"!="" {
				replace _weightvar = `weightvar'
			}
			
			*Normalize weights so their sum equals N
			sum _weightvar, meanonly
			replace _weightvar = _weightvar * (`N' / r(sum))
			
			*Get the independent variables for the X matrix
			local names: colnames e(b)
			local n = 1
			foreach name of local names {
				cap drop _tmp`n'
				gen _tmp`n' = `name'
				sum _tmp`n', meanonly
				if r(mean) != 0 {
					if `n' == 1 local variables "_tmp`n'"
					else        local variables "`variables' _tmp`n'"
				}
				else {
					drop _tmp`n'
				}
				local ++n
			}
			
			*Store matrices for Python access
			matrix _bds = `tmp_w_mat'
			matrix _min_pval = `min_pval'
			matrix _max_pval = `max_pval'
			
			*Clear Python environment and any existing result matrices
			python clear
			cap matrix drop _orig _costs _bds _min_pval _max_pval
			noi python script "`c(sysdir_plus)'py/p_value_cost_minimizer.py"
			
			*Get results from Python and store in temp matrices
			tempname costs_pval_python orig_pval_python
			matrix `costs_pval_python' = _costs
			matrix `orig_pval_python' = _orig
			
			*Clean up temporary matrices from Python
			matrix drop _costs _orig _bds _min_pval _max_pval
			
			restore
		}
	}
	
	*=====================================
	*4. Loop over levels of c to find reversal conditions (only if pythonno specified)
	*=====================================
	
	if "`pythonno'" != "" {
		
		*-------------------------------------
		*4.1 Initialising things 
		*-------------------------------------

		*Tell the user how many regressions need to be run
		local N_c = ((`end' - `start')/`precision')+1
		noisily dis ""
		noisily dis "Total number of c-values to check: `N_c'. Progress:"
		
		*A counter to keep track of things
		local iter = 1
		
		*-------------------------------------
		*4.2 Begin the loop 
		*-------------------------------------

		forvalues c=`start'(`precision')`end'{
			
			*-------------------------------------
			*4.3 Get minus sign before the exp() function for negative c
			*-------------------------------------

			if `c' < 0 	local minus "-"
			if `c' > 0	local minus ""
			
			*=====================================
			*4.4 Run the fast routine
			*=====================================

			if "`fast'"=="fast" {
				
				*-------------------------------------
				*4.4.1 Get a rescaled version of the original variable
				*-------------------------------------

				tempvar trans_depvar
				if (`c' < 0.0000001 & `c' > -0.0000001) gen `trans_depvar' = `depvar' // I.e. if c==0
				else {
					gen `trans_depvar' = ((`minus'exp(`depvar'*`c') - `minus'exp(`scale_min'*`c')) / (`minus'exp(`scale_max'*`c') - `minus'exp(`scale_min'*`c')))*(`scale_max'-`scale_min') + `scale_min'
				}

				*-------------------------------------
				*4.4.2 Get a matrix to record differences in labels
				*-------------------------------------

				matrix `dlabels' = J(`nrows_d_result',`nrows_d_result',0)
				local n = 0
				levelsof `trans_depvar', local(levels_trans_depvar)
				foreach l of local levels_trans_depvar {
					if `n'==0 {
						local k = `l' 	// k here denotes the previous label. Do nothing in the first iteration and just update k for the next iteration.
					}
					else { 				// produces a diagonal matrix with the right dimensions to multiply.
						capture matrix `dlabels'[`n',`n'] = `k'-`l'  									// Again weird condition because of rounding errors from floating points.
						local k = `l' 	// Store l in k for the next iteration. 
					}
					local ++n
				}
				
				*-------------------------------------
				*4.4.3 Use the results from regressions of hd to arrive at the transformed coefficients. 
				*-------------------------------------

				*Take the dmat result and multiply each element by the matrix recording differences of labels
				matrix `tmp1' = `dlabels' * `tmp_w_mat'
				
				*Take the column-wise sum:
				matrix `tmp1' = J(1,rowsof(`tmp1'),1)*`tmp1' 
				
				*Store result in the b_result matrix
				if `c'==`start' matrix `b_result' = (`tmp1')
				if `c'!=`start' matrix `b_result' = (`b_result' \ `tmp1')
				
			}
				
			*=====================================
			*4.5 Run the slow routine (only needed when interest is in p-values.)
			*=====================================

			if "`fast'" == "" {
				
				*-------------------------------------
				*4.5.1 Create transformed depvar and run the regression
				*-------------------------------------
				
				*Create transformed depvar			
				tempvar trans_depvar
				if (`c' < 0.0000001 & `c' > -0.0000001) gen `trans_depvar' = `depvar' 
				else {
					gen double `trans_depvar' = ((`minus'exp(`depvar'*`c') - `minus'exp(`scale_min'*`c')) / (`minus'exp(`scale_max'*`c') - `minus'exp(`scale_min'*`c')))*(`scale_max'-`scale_min') + `scale_min'
				}
				*Get the command to run
				local to_run = subinword("`full_command'", "`depvar'", "`trans_depvar'", 1) 	// changes the dependent variable of the command originally run.
				
				*Run the regression
				`to_run'

				*Save results
				matrix `tmp1' = e(b)
				matrix `tmp2' = vecdiag(e(V))	
							
				*-------------------------------------
				*4.5.2 Get p-values and put coefficients into matrices  
				*-------------------------------------
				
				*Put coefficients into a results-matrix		
				if `c'==`start' matrix `b_result' = (`tmp1')			
				if `c'!=`start' matrix `b_result' = (`b_result' \ `tmp1')

				*Get p-vals and put into a results-matrix
				tempname tmp_p
				
				if `ncluster'!=. capture mata : st_matrix("`tmp_p'", 2*(J(1,`P',1)-t((`ncluster'-1), abs(st_matrix("`tmp1'") :/ sqrt(st_matrix("`tmp2'")))))) 			
				if `ncluster'==. capture mata : st_matrix("`tmp_p'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`tmp1'") :/ sqrt(st_matrix("`tmp2'")))))) 

				if `c'!=`start' matrix `p_result' = (`p_result' \ `tmp_p')
				if `c'==`start' matrix `p_result' = (`tmp_p')

				
			}
			
			*Display progress
			noisily _dots `iter' 0
			local ++iter
		}
			
		*=====================================
		*4.6 Use transformed estimates from any of the routines to find points of reversal across c
		*=====================================
		
		preserve
		clear
		
		svmat double `b_result', names(variable)
		qui gen tmpvar = _n
		qui tsset tmpvar
		foreach var of varlist variable1-variable`P' {
			if "`var'" == "tmpvar" continue, break  
			qui replace `var' = `var'>`revpoint' 		// 1 if positive, 0 otherwise
			qui gen double new`var' = d.`var' 	// missing in the first instance, then -1 if pos to neg, and 1 if neg to pos, 0 if the same.
		}
		keep newvariable1-newvariable`P'
		tempname r_result
		mkmat newvariable1-newvariable`P', matrix(`r_result')
		
		*=====================================
		*4.7 Use p-value estimates from to find points of significant reversal across c (only for slow routine)
		*=====================================
		
		if "`fast'"=="" {
			clear
			svmat double `p_result', names(variable)	
			qui gen tmpvar = _n
			qui tsset tmpvar	
			foreach var of varlist variable1-variable`P' {
				if "`var'" == "tmpvar" continue, break  
				qui replace `var' = (`var'-`critval')<0 	// 1 if significant, 0 if insifnificant
				qui gen double new`var' = d.`var' 			// missing in the first instance, then -1 if significant to insignificant, and 1 if insignificant to significant, 0 if the same.
			}		
			keep newvariable1-newvariable`P'
			tempname y_result
			mkmat newvariable1-newvariable`P', matrix(`y_result')
		}
		restore
	
	} // End of pythonno conditional
	
	*=====================================
	*5. Clean-up matrices
	*=====================================

	*-------------------------------------
	*5.1 Get the original variable names as column names 
	*-------------------------------------

	* Only set column names for matrices that exist
	if "`pythonno'" == "" {
		matrix colnames `label_mat'=`explanatory_vars'
		matrix colnames `cost_mat'=`explanatory_vars'
	}
	
	if "`pythonno'" != "" {
		matrix colnames `b_result'=`:colnames `tmp_w_mat'' 
		matrix colnames `r_result'=`:colnames `tmp_w_mat'' 
		if "`fast'"=="" matrix colnames `p_result'=`:colnames `tmp_w_mat'' 
		if "`fast'"=="" matrix colnames `y_result'=`:colnames `tmp_w_mat''
	}
	
	matrix colnames `d_result'=`:colnames `tmp_w_mat'' 
	matrix colnames `hd_p_vals'=`:colnames `tmp_w_mat'' 
	
	*-------------------------------------
	*5.2 Attach the corresponding value of c to the results matrices (only if pythonno)
	*-------------------------------------
	
	if "`pythonno'" != "" {
		local n = 1
		tempname c_vals
		matrix `c_vals' = J(rowsof(`b_result'),1,.)
		forvalues c=`start'(`precision')`end'{
			capture matrix `c_vals'[`n',1] = `c' // not sure why there's a capture here. 
			local ++n
		}
		
		matrix `b_result' = (`b_result', `c_vals')
		matrix `r_result' = (`r_result', `c_vals') 	
		
		if "`fast'"=="" matrix `p_result' = (`p_result', `c_vals') 
		if "`fast'"=="" matrix `y_result' = (`y_result', `c_vals')
	} 

	*=====================================
	*6. Preparare for main display to user
	*=====================================
	
	*-------------------------------------
	*6.1 Get minimum c-value (only if pythonno option specified)
	*-------------------------------------

	if "`pythonno'" != "" {
		preserve
		clear
		local c_num = colsof(`r_result')
		local c_num1 = `c_num' - 1
		local c_num2 = `c_num' - 2
		
		svmat `r_result', names(tmp)
		gen orig_c = tmp`c_num' 
		
		forvalues j=1/`c_num' {
			replace tmp`j' = abs(tmp`j') 	
		}
		forvalues j=1/`c_num1' {
			sum tmp`c_num' if tmp`j'==1
			sum orig_c if tmp`c_num'==r(min) & tmp`j'==1
			replace tmp`j' = r(mean) 	
		}
		keep if _n==1
		tempname min_c_value
		mkmat tmp1-tmp`c_num1', matrix(`min_c_value')
		restore
	}
	
	if "`pythonno'" != "" & "`pvalue'" != "" & "`fast'" == "" {

		*-------------------------------------
		*6.1.1 Get minimum c-value to render significant coefficient insignificant
		*-------------------------------------

		/*
		Gameplan:
		Step 1: Only keep the matrix for positive c. Find first -1.
		Step 2: Only keep the matrix for negative c. Take absolute value of c, sort on c. Find first 1.
		Step 3: Combine and find smallest c. 
		*/ 
		
		*Step 1:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' >= 0

		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==-1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp_value1
		mkmat tmp1-tmp`c_num2', matrix(`min_cp_value1')
		restore
		
		
		*Step 2:
		if `start' < 0 {
			preserve
			clear
			svmat `y_result', names(tmp)
			keep if tmp`c_num' < 0
			replace tmp`c_num' = abs(tmp`c_num')
			sort tmp`c_num'
			forvalues j=1/`c_num2' {
				sum tmp`c_num' if tmp`j'==1
				replace tmp`j' = r(min) 	
			}
			keep if _n==1
			tempname min_cp_value2
			mkmat tmp1-tmp`c_num2', matrix(`min_cp_value2')
			restore
		}
		
		*Step 3:
		tempname combined
		if `start' < 0 	matrix `combined' = `min_cp_value1' \ `min_cp_value2'
		else 			matrix `combined' = `min_cp_value1'
		preserve 
		clear
		svmat `combined', names(tmp)
		gen sign = 1 in 1
		if `start' < 0 replace sign = -1 in 2
		
		forvalues j=1/`c_num2' {
			sum tmp`j' 
			sum sign if tmp`j' == r(min)
			local sign = r(mean)
			sum tmp`j' 		
			replace tmp`j' = r(min) * `sign'	
		}
		keep if _n==1
		tempname min_cp_value
		mkmat tmp1-tmp`c_num2', matrix(`min_cp_value')
		restore
		
		*-------------------------------------
		*6.1.2 Get minimum c-value to render insignificant coefficient significant
		*-------------------------------------
		
		/*
		Gameplan:
		Step 1: Only keep the matrix for positive c. Find first 1.
		Step 2: Only keep the matrix for negative c. Take absolute value of c, sort on c. Find first -1.
		Step 3: Combine and find smallest c. 
		*/ 

		*Step 1:
		preserve
		clear
		svmat `y_result', names(tmp)
		keep if tmp`c_num' >= 0

		forvalues j=1/`c_num2' {
			sum tmp`c_num' if tmp`j'==1
			replace tmp`j' = r(min) 	
		}
		keep if _n==1
		tempname min_cp2_value1
		mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value1')
		restore
		
		
		*Step 2:
		if `start' < 0 {
			preserve
			clear
			svmat `y_result', names(tmp)
			keep if tmp`c_num' < 0
			replace tmp`c_num' = abs(tmp`c_num')
			sort tmp`c_num'
			forvalues j=1/`c_num2' {
				sum tmp`c_num' if tmp`j'==-1
				replace tmp`j' = r(min) 	
			}
			keep if _n==1
			tempname min_cp2_value2
			mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value2')
			restore
		}
		
		*Step 3:
		tempname combined
		if `start' < 0 	matrix `combined' = `min_cp2_value1' \ `min_cp2_value2'
		else			matrix `combined' = `min_cp2_value1' 
		preserve 
		clear
		svmat `combined', names(tmp)
		gen sign = 1 in 1
		if `start' < 0 replace sign = -1 in 2
		
		forvalues j=1/`c_num2' {
			sum tmp`j' 
			sum sign if tmp`j' == r(min)
			local sign = r(mean)
			sum tmp`j' 		
			replace tmp`j' = r(min) * `sign'	
		}
		keep if _n==1
		tempname min_cp2_value
		mkmat tmp1-tmp`c_num2', matrix(`min_cp2_value')
		restore	
	}
	*-------------------------------------
	*6.2 Get matrices into the right format, find original p-values, construct display matrix, implement the "keep()" option.
	*-------------------------------------
	
	**Get the right column number (removing the constant)		
	local npcols = colsof(`orig_b') - 1 
	
	**Cut to right number of columns
	*Original coefficients.
	tempname orig_b_display
	matrix `orig_b_display' = `orig_b'[1...,1..`npcols']
	
	*Minimum and maximum p-value (only if computed)
	if "`pvalue'" != "" {
		matrix `max_pval' = `max_pval'[1...,1..`npcols'] 
		matrix `min_pval' = `min_pval'[1...,1..`npcols'] 
	}
	
	*Minimum c-value (only if computed)
	if "`pythonno'" != "" {
		matrix `min_c_value' = `min_c_value'[1...,1..`npcols']
	}

	*Original p-value (only compute if pvalue option specified)
	if "`pvalue'" != "" {
		tempname orig_ses
		matrix `orig_ses' = vecdiag(`orig_v')
		tempname orig_p_display
		if `ncluster'!=. mata : st_matrix("`orig_p_display'", 2*(J(1,`P',1)-t((`ncluster'-1), abs(st_matrix("`orig_b'") :/ sqrt(st_matrix("`orig_ses'")))))) 
		if `ncluster'==. mata : st_matrix("`orig_p_display'", 2*(J(1,`P',1)-t((`N'-`P'), abs(st_matrix("`orig_b'") :/ sqrt(st_matrix("`orig_ses'")))))) 	
		matrix `orig_p_display' = `orig_p_display'[1...,1..`npcols']
	} 

	*-------------------------------------
	*6.2.1 Calculate upper and lower bounds for each coefficient from bd matrix 
	*-------------------------------------
	
	**Calculate bounds from bd coefficients (stored in tmp_w_mat)
	tempname coeff_upper_bounds coeff_lower_bounds
	
	*Get the right column number (removing the constant)
	local bd_npcols = colsof(`tmp_w_mat') - 1
	
	*Cut to right number of columns for bd matrix 
	tempname tmp_w_mat_display
	matrix `tmp_w_mat_display' = `tmp_w_mat'[1...,1..`bd_npcols']
	
	*Calculate upper bounds: -1 * minimum hd coefficient * scale range
	preserve
	clear
	svmat `tmp_w_mat_display', names("bd")
	ds
	local bd_vars = "`r(varlist)'"
	foreach var in `bd_vars' {
		sum `var'
		replace `var' = -r(min) * (`scale_max' - `scale_min')
	}
	keep if _n == 1
	mkmat `bd_vars', matrix(`coeff_upper_bounds')
	restore
	
	*Calculate lower bounds: -1 * maximum hd coefficient * scale range
	preserve  
	clear
	svmat `tmp_w_mat_display', names("bd")
	ds
	local bd_vars = "`r(varlist)'"
	foreach var in `bd_vars' {
		sum `var'
		replace `var' = -r(max) * (`scale_max' - `scale_min')
	}
	keep if _n == 1
	mkmat `bd_vars', matrix(`coeff_lower_bounds')
	restore
	

	*-------------------------------------
	*6.3 Assemble the matricies to be displayed
	*-------------------------------------
	
	tempname _to_display
	matrix `_to_display' = `orig_b_display'
	matrix `_to_display' = `_to_display' \ `coeff_lower_bounds'
	matrix `_to_display' = `_to_display' \ `coeff_upper_bounds'
		
	* Default behavior: Show only Python cost (if it exists)
	if "`pythonno'"=="" {
		capture confirm matrix `cost_mat'
		if _rc == 0 {
			matrix `_to_display' = `_to_display' \ `cost_mat'
		}
	}
	
	* pythonno behavior: Show minimum c-value (if it exists)
	if "`pythonno'" != "" {
		capture confirm matrix `min_c_value'
		if _rc == 0 {
			matrix `_to_display' = `_to_display' \ `min_c_value'
		}
	}
	
	* pvalue option: Show p-value statistics
	if "`pvalue'" != "" {
		* Show original p-value, then p-value bounds (from hd regressions)
		matrix `_to_display' = `_to_display' \ `orig_p_display' \ `min_pval' \ `max_pval'
	}
	
	* Show cost for p-value reversals (python case)
	if "`pvalue'" != "" &  "`pythonno'" == "" {
		mat `costs_pval_python' = `costs_pval_python''
		local ncols_for_pvals = colsof(`costs_pval_python') - 1
		mat `costs_pval_python' = `costs_pval_python'[1,1..`ncols_for_pvals']	
		matrix `_to_display' = `_to_display' \ `costs_pval_python'
		
		* Show original p-values from Python (for debugging)
		// mat `orig_pval_python' = `orig_pval_python''
		// mat `orig_pval_python' = `orig_pval_python'[1,1..`ncols_for_pvals']
		// matrix `_to_display' = `_to_display' \ `orig_pval_python'
	}
	
	* Show cost for p-value reversals (non-python case)
	if "`pvalue'" != "" &  "`pythonno'" != "" {
		matrix `_to_display' = `_to_display' \ `min_cp_value'
	}
	
	*Implement keep option
	if "`keep'" != "" matselrc `_to_display' `_to_display', c(`keep') 
	
	*-------------------------------------
	*6.3.1 Get variable names for r() returns after keep() operation
	*-------------------------------------
	
	* Get variable names after keep() has been applied
	local var_names : colnames `_to_display'
	
	*-------------------------------------
	*6.4 Get the labels and notes
	*-------------------------------------
	
	*Labels  
	local display_labels "Coef"
	local display_labels "`display_labels'" "Min.coef" "Max.coef"  
	
	* Default behavior: Show only Python cost (if it exists)
	if "`pythonno'"=="" {
		capture confirm matrix `cost_mat'
		if _rc == 0 {
			local display_labels "`display_labels'" "Min.cost"
		}
	}
	
	* pythonno behavior: Show minimum c-value (if it exists)
	if "`pythonno'" != "" {
		capture confirm matrix `min_c_value'
		if _rc == 0 {
			local display_labels "`display_labels'" "Min.c"
		}
	}
	
	* pvalue option: Show p-value statistics
	if "`pvalue'" != "" {
		local display_labels "`display_labels'" "P-val" "Min.p-val" "Max.p-val"
		
		* Add cost label for p-value target (only if computed by Python for non-clustered SEs)
		if "`clustvar'" == "" {
			capture confirm matrix `costs_pval_python'
			if _rc == 0 {
				local display_labels "`display_labels'" "Min.cost.sig."
			}
			// capture confirm matrix `orig_pval_python'
			// if _rc == 0 {
			//	local display_labels "`display_labels'" "Py.Pval"
			// }
		}
	}
	
	* Show minimum c-value for significance change label (displayed last, if it exists)
	if "`pythonno'" != "" {
		capture confirm matrix `min_cp_value'
		if _rc == 0 {
			local display_labels "`display_labels'" "Min.c.sig."
		}
	}
	
	*Notes
	local notes "Note: {bf:Missing values imply one of the following:} `=char(13)' (1) That no original coefficient was estimated. `=char(13)' (2) That no reversal is possible. `=char(13)' (3) That the minimum reversing c-value falls outside the search range. `=char(13)'"
	if "`pythonno'" != "" & "`pvalue'" != "" & "`fast'" == "" {
		local notes `notes' "(4) That coefficients cannot be made significant or insignificant within the search range."
	}
	
	
	*=====================================
	*7. Display to user 
	*====================================
	
	noi {
		
		*-------------------------------------
		*7.1 Main display
		*-------------------------------------
	
		dis ""
		dis ""		
		dis "{bf:Results:}"
		esttab matrix(`_to_display', fmt(3) transpose), mtitles("") ///
		collabels("`display_labels'") ///
		modelwidth(13) note( `notes') //
		
		*-------------------------------------
		*7.2 Display column explanations
		*-------------------------------------
		
		dis ""
		dis "{bf:Column explanations:}"
		dis "Coef: Original coefficient from fitted model"		
		dis "Min.coef: Lower bound of coefficient across hd transformations"
		dis "Max.coef: Upper bound of coefficient across hd transformations"
		
		if "`pythonno'"=="" {
			capture confirm matrix `cost_mat'
			if _rc == 0 {
				dis "Min.cost: Minimum cost for coefficient sign reversal (Python optimization)"
			}
		}
		
		if "`pythonno'" != "" {
			capture confirm matrix `min_c_value'
			if _rc == 0 {
				dis "Min.c: Minimum c-value for coefficient sign reversal (exponential transformation)"
			}
		}
		
		if "`pvalue'" != "" {
			dis "P-val: Original p-value from fitted model"
			dis "Min.p-val: Minimum p-value across hd transformations"
			dis "Max.p-val: Maximum p-value across hd transformations"
			capture confirm matrix `costs_pval_python'
			if _rc == 0 {
					dis "Min.cost.sig.: Minimum cost for statistical significance reversal"
			}
		}
		
		if "`pythonno'" != "" {
			capture confirm matrix `min_cp_value'
			if _rc == 0 {
				dis "Min.c.sig.: Minimum c-value for statistical significance reversal"
			}
		}
	}


	*=====================================
	*8. Return results in r()
	*=====================================
	
	*-------------------------------------
	*8.1 Return main result matrix (complete displayed table)
	*-------------------------------------
	
	tempname return_result
	matrix `return_result' = `_to_display'
	matrix colnames `return_result' = `var_names'
	return matrix result `return_result'
	
	*-------------------------------------
	*8.3 Return coefficient matrices
	*-------------------------------------
	
	* r(b) - Original coefficients
	tempname return_b
	matrix `return_b' = `orig_b_display'
	matrix colnames `return_b' = `var_names'
	return matrix b `return_b'
	
	* r(minb) and r(maxb) - Coefficient bounds (only if bounds computed)
	tempname return_minb return_maxb
	matrix `return_minb' = `coeff_lower_bounds'
	matrix `return_maxb' = `coeff_upper_bounds'
	matrix colnames `return_minb' = `var_names'
	matrix colnames `return_maxb' = `var_names'
	return matrix minb `return_minb'
	return matrix maxb `return_maxb'
	
	
	*-------------------------------------
	*8.4 Return minimum reversing c-values
	*-------------------------------------
	
	* r(minc) - Minimum c-values for coefficient reversal
	if "`pythonno'" != "" {
		capture confirm matrix `min_c_value'
		if _rc == 0 {
			tempname return_minc
			matrix `return_minc' = `min_c_value'
			matrix colnames `return_minc' = `var_names'
			return matrix minc `return_minc'
		}
	}
	
	* Return Python cost matrix if available
	if "`pythonno'" == "" {
		capture confirm matrix `cost_mat'
		if _rc == 0 {
			tempname return_cost
			matrix `return_cost' = `cost_mat'
			matrix colnames `return_cost' = `var_names'
			return matrix cost `return_cost'
		}		
	}
	
	*-------------------------------------
	*8.5 Return p-value matrices
	*-------------------------------------
	
	if "`pvalue'" != "" {
		* r(p) - Original p-values
		tempname return_p
		matrix `return_p' = `orig_p_display'
		matrix colnames `return_p' = `var_names'
		return matrix p `return_p'
		
		* r(minp) and r(maxp) - P-value bounds (from hd regressions)
		tempname return_minp return_maxp
		matrix `return_minp' = `min_pval'
		matrix `return_maxp' = `max_pval'
		matrix colnames `return_minp' = `var_names'
		matrix colnames `return_maxp' = `var_names'
		return matrix minp `return_minp'
		return matrix maxp `return_maxp'
		
		* r(costp) - Cost for achieving target p-value (only if computed by Python for non-clustered SEs)
		if "`clustvar'" == "" {
			capture confirm matrix `costs_pval_python'
			if _rc == 0 {
				tempname return_costp
				matrix `return_costp' = `costs_pval_python'
				matrix colnames `return_costp' = `var_names'
				return matrix costp `return_costp'
			}
		}
	}
	
	*-------------------------------------
	*8.6 Return minimum c-values for significance reversal
	*-------------------------------------
	
	* r(mincp) - Minimum c-values for significance change
	if "`pythonno'" != "" {
		capture confirm matrix `min_cp_value'
		if _rc == 0 {
			tempname return_mincp
			matrix `return_mincp' = `min_cp_value'
			matrix colnames `return_mincp' = `var_names'
			return matrix mincp `return_mincp'
		}
	}
	
	*-------------------------------------
	*8.7 Return additional matrices for advanced users
	*-------------------------------------
	
	* Return internal matrices
	tempname return_d return_hdp
	matrix `return_d' = `d_result'
	matrix `return_hdp' = `hd_p_vals'
	matrix colnames `return_d' = `var_names'
	matrix colnames `return_hdp' = `var_names'
	return matrix d `return_d' 
	return matrix hdp `return_hdp'
	
	if "`pythonno'" != "" {
		capture confirm matrix `b_result'
		if _rc == 0 {
			return matrix b_full `b_result' 
			return matrix r_full `r_result'
		}
		
		* Only return p-value matrices if slow routine was used
		if "`pvalue'" != "" & "`fast'" == "" {
			capture confirm matrix `p_result'
			if _rc == 0 {
				return matrix p_full `p_result' 
				return matrix y_full `y_result'
			}
		}
	} 
		
	*=====================================
	*9. Restore the original model.
	*=====================================
	
	cap drop _hd_residual*
	matrix drop _labels_depvar
	qui estimates restore `prevmodel'
	
	}	// ends the qui condition
	
end
