{smcl}
help for {hi:coeff_reverser} version 0.2 {right: (Caspar Kaiser & Anthony Lepinteur)}
{hline}
{title:Coefficient Reversal Analysis for Ordered Response Data}

{title:Syntax}

{p 4 4 2}
{cmd:coeff_reverser}
[{cmd:,} {it:options}]

{synoptset tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Stata-based search {help coeff_reverser##opt_stata:[+]}}
{synopt:{cmd:max_attempts(}{it:integer}{cmd:)}}Cap on the number of faces evaluated by the analytical active-set search (default: 50000){p_end}

{syntab:Python-based search {help coeff_reverser##opt_python:[+]}}
{synopt:{opt python}}Use Python numerical optimisation instead of the default analytical active-set search{p_end}
{synopt:{opt pvalue}}Display p-value statistics (min/max p-values, and costs for significance changes). Under development.{p_end}
{synopt:{cmd:revpoint(}{it:real}{cmd:)}}Specifies the target value for sign reversal (default: 0){p_end}
{synopt:{cmd:critval(}{it:real}{cmd:)}}Specifies the level for statistical significance (default: 0.05){p_end}
{synopt:{cmd:alpha(}{it:real}{cmd:)}}Specifies the alpha parameter for the cost function (default: 2){p_end}
{synopt:{opt theil}}Use normalized Theil index as cost function{p_end}

{syntab:Output options {help coeff_reverser##opt_output:[+]}}
{synopt:{cmd:keep(}{it:string}{cmd:)}}Specifies list of variables to keep in displayed results table{p_end}

{synoptline}
{p 4 4} {cmd:coeff_reverser} requires that a linear regression has been run and is stored in {cmd:e()}. {p_end}
{p 4 4} When {opt python} is specified, the command requires Python to be available, and SciPy and NumPy to be available in Stata's Python environment. {p_end}
{p 4 4} Also see {help mrs_reverser}. {p_end}


{marker introduction}{...}
{title:Introduction}

{p 4 4}{cmd:coeff_reverser} analyzes whether coefficient signs can be reversed under some positive monotonic transformations of the dependent variable.
This implements the methodology developed in Kaiser & Lepinteur (2025) for testing the robustness of empirical results to departures from the linearity assumption in survey scales.

{p 4 4}See Kaiser & Lepinteur (2025; {browse "https://arxiv.org/abs/2507.16440v1":https://arxiv.org/abs/2507.16440}) for details.{p_end}

{marker description}{...}
{title:Description}

{p 4 4}{cmd:coeff_reverser} takes a linear regression and determines whether coefficient signs can be reversed through monotonic transformations of the dependent variable.

The command provides several outputs:

{p 6 8}• Whether each coefficient's sign can be reversed{p_end}
{p 6 8}• Minimum 'cost' (i.e. departure from linearity) required for reversal{p_end}
{p 6 8}• Upper and lower bounds for coefficients under any transformations{p_end}
{p 6 8}• P-value bounds showing how statistical significance can change ({opt python} only){p_end}
{p 6 8}• Minimum 'cost' required to achieve a certain p-value ({opt python} only){p_end}


Two computational approaches are available:

{p 6 8}• {bf:Analytical active-set search} (default): Uses an analytical active-set search written in Stata/Mata to find 'least non-linear' transformations needed to achieve sign reversal under the SD-based cost function. {p_end}
{p 6 8}• {bf:Python optimisation}: Uses Python to numerically minimise a cost function as described in Kaiser & Lepinteur (2025). Required when working with cost functions other than the SD-based default (e.g. {opt theil}), with non-zero reversal targets ({opt revpoint}), or for p-value bounds ({opt pvalue}). {p_end}

{title:Options}
{marker opt_stata}{...}
{dlgtab:Stata-based search}

{p 4 4} {cmd:max_attempts(}{it:integer}{cmd:)} caps the number of faces evaluated by the analytical active-set search. Default is 50000. If the cap is reached before the search completes, the command reports the best feasible cost found so far and issues a warning. Ignored when {opt python} is specified.

{marker opt_python}{...}
{dlgtab:Python-based search}

{p 4 4} {opt python} switches from the default analytical active-set search to Python-based numerical optimisation. Required to use {opt pvalue}, {opt revpoint}, {opt critval}, {opt alpha}, or {opt theil}.

{p 4 4} {opt pvalue} displays p-value statistics including original p-values, minimum and maximum p-values achievable through transformations, and minimum costs needed to change statistical significance. Only relevant when {opt python} is specified. Currently only works after running {cmd:reg} and only for standard and 'robust' standard errors (and {bf:not} for clustered SEs). Also does {bf:not} work with factor variables.

{p 4 4} {cmd:revpoint(}{it:real}{cmd:)} specifies the target value for coefficient reversal. Default is 0 (sign reversal).
For example, {cmd:revpoint(0.5)} checks if coefficients can be transformed to equal 0.5, and, if so, at what cost. Only relevant when {opt python} is specified.

{p 4 4} {cmd:critval(}{it:real}{cmd:)} sets the significance level for statistical tests. Default is 0.05. Only relevant when {opt pvalue} is specified.

{p 4 4} {cmd:alpha(}{it:real}{cmd:)} controls the cost function sensitivity. Higher values penalize uneven transformations more heavily.
The cost function is: cost = (var/maxvar)^(1/alpha). Default is alpha=2. Only relevant when {opt python} is specified.

{p 4 4} {opt theil} uses the normalized Theil inequality index as the cost function instead of the alpha-based variance cost function. Only relevant when {opt python} is specified.

{marker opt_output}{...}
{dlgtab:Output options}

{p 4 4} {cmd:keep(}{it:string}{cmd:)} restricts the output to specific variables from the original regression.

{marker examples}{...}
{title:Examples}

{p 4 4}Basic reversal analysis after linear regression:{p_end}
{p 8 12}{inp:. sysuse auto}{p_end}
{p 8 12}{inp:. regress rep78 price mpg headroom displacement gear_ratio}{p_end}
{p 8 12}{inp:. coeff_reverser}{p_end}

{p 4 4}Tighten or loosen the cap on the analytical search:{p_end}
{p 8 12}{inp:. coeff_reverser, max_attempts(100000)}{p_end}

{p 4 4}Switch to Python optimisation:{p_end}
{p 8 12}{inp:. coeff_reverser, python}{p_end}

{p 4 4}Include p-value bounds (requires {opt python}):{p_end}
{p 8 12}{inp:. coeff_reverser, python pvalue}{p_end}

{p 4 4}Use Python optimisation with Theil cost function:{p_end}
{p 8 12}{inp:. coeff_reverser, python theil}{p_end}

{p 4 4}Check reversal to specific target value (here =0.2):{p_end}
{p 8 12}{inp:. coeff_reverser, python revpoint(0.2)}{p_end}

{p 4 4}Only display specific variables:{p_end}
{p 8 12}{inp:. coeff_reverser, keep(income education)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:coeff_reverser} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Main results}{p_end}
{synopt:{cmd:r(result)}}complete displayed table with all computed statistics{p_end}

{p2col 5 20 24 2: Coefficient matrices}{p_end}
{synopt:{cmd:r(b)}}original coefficient estimates{p_end}
{synopt:{cmd:r(minb)}}lower bounds for coefficients{p_end}
{synopt:{cmd:r(maxb)}}upper bounds for coefficients{p_end}

{p2col 5 20 24 2: Cost matrices}{p_end}
{synopt:{cmd:r(cost)}}minimum transformation cost for coefficient sign reversal{p_end}

{p2col 5 20 24 2: Analytical search diagnostics (default mode only)}{p_end}
{synopt:{cmd:r(Delta_star)}}optimal gap vectors achieving the minimum-cost reversal{p_end}
{synopt:{cmd:r(A_bits)}}bitmask encoding the active set at the optimum{p_end}
{synopt:{cmd:r(n_attempts)}}number of faces evaluated by the search{p_end}
{synopt:{cmd:r(hit_limit)}}indicator (per coefficient) for whether the search hit {opt max_attempts}{p_end}

{p2col 5 20 24 2: P-value matrices (with {opt python} and {opt pvalue})}{p_end}
{synopt:{cmd:r(p)}}original p-values from fitted model{p_end}
{synopt:{cmd:r(minp)}}minimum p-values across transformations{p_end}
{synopt:{cmd:r(maxp)}}maximum p-values across transformations{p_end}
{synopt:{cmd:r(costp)}}transformation costs for significance reversal{p_end}

{p2col 5 20 24 2: Advanced matrices}{p_end}
{synopt:{cmd:r(d)}}reversal indicators for each coefficient{p_end}
{synopt:{cmd:r(hdp)}}p-values from hd transformations{p_end}

{title:References}

{p 4 4} Kaiser, Caspar and Anthony Lepinteur. (2025). Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data. Available at {browse "https://arxiv.org/abs/2507.16440v1":https://arxiv.org/abs/2507.16440}.

{title:Author}

{p 4 4}  Caspar Kaiser {p_end}
{p 4 4}  Warwick Business School, University of Warwick {p_end}
{p 4 4}  caspar.kaiser@wbs.ac.uk {p_end}

{p 4 4}  Anthony Lepinteur {p_end}
{p 4 4}  University of Luxembourg {p_end}
{p 4 4}  anthony.lepinteur@uni.lu {p_end}

{p 4 4} {cmd:coeff_reverser} is under active development. Please report bugs or feature requests to the authors.
