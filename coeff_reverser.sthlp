{smcl}
help for {hi:coeff_reverser} version 0.1 {right: (Caspar Kaiser & Anthony Lepinteur)}
{hline}
{title:Coefficient Reversal Analysis for Ordered Response Data}

{title:Syntax}

{p 4 4 2}
{cmd:coeff_reverser}
[{cmd:,} {it:options}]

{synoptset tabbed}{...}
{synopthdr}
{synoptline}
{syntab:General {help coeff_reverser##opt_general:[+]}}
{synopt:{opt pythonno}}For use when python is not available. Searches over exponential transformations instead of using the cost-function approach. {p_end}
{synopt:{opt pvalue}}Display p-value statistics (min/max p-values, and costs for significance changes){p_end}
{synopt:{cmd:revpoint(}{it:real}{cmd:)}}Specifies the target value for sign reversal (default: 0){p_end}
{synopt:{cmd:critval(}{it:real}{cmd:)}}Specifies the level for statistical significance (default: 0.05){p_end}

{syntab:Cost-function options {help coeff_reverser##opt_search:[+]}}
{synopt:{cmd:alpha(}{it:real}{cmd:)}}Specifies the alpha parameter for the cost function (default: 2){p_end}
{synopt:{opt theil}}Use normalized Theil index as cost function (overrides {cmd:alpha} option){p_end}

{syntab:Exponential function search options (applies when specifying {cmd:pythonno}) {help coeff_reverser##opt_search:[+]}}
{synopt:{cmd:start(}{it:real}{cmd:)}}Smallest value of c over which to search (default: -2){p_end}
{synopt:{cmd:end(}{it:real}{cmd:)}}Largest value of c over which to search (default: 2){p_end}
{synopt:{cmd:precision(}{it:real}{cmd:)}}Grid precision for c values (default: 0.1){p_end}

{syntab:Output options {help coeff_reverser##opt_output:[+]}}
{synopt:{cmd:keep(}{it:string}{cmd:)}}Specifies list of variables to keep in displayed results table{p_end}

{synoptline}
{p 4 4} {cmd:coeff_reverser} requires that a regression has been run and is stored in {cmd:e()}. {p_end}
{p 4 4} Unless pythonno is specified, this command requires Python to be available, and SciPy and NumPy to be available in Stata's Python environment. {p_end}

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
{p 6 8}• P-value bounds showing how statistical significance can change{p_end}
{p 6 8}• Minimum 'cost' required to achieve a certain p-value{p_end}


Two computational approaches are available:

{p 6 8}• {bf:Python optimization}: Uses Python to find transformations minimizing a cost function as described in Kaiser & Lepinteur (2025). This finds 'least non-linear' transformations needed to achieve sign reversal. {p_end}
{p 6 8}• {bf:Exponential transformation search}: Should Python be unavailable, searches over exponential transformations of the form f(depvar)=exp(depvar*c) (for positive c) or f(depvar)=-exp(depvar*c) (for negative c). This follows the approach of Bond & Lang (2019) and Kaiser & Vendrik (2023). This approach is much more restrictive and only finds least non-linear reversals within the exponential class.{p_end}

{title:Options}
{marker opt_general}{...}
{dlgtab:General}

{p 4 4} {opt pythonno} for use when Python is not available. Searches over exponential transformations of the form f(depvar)=exp(depvar*c) (if c>0) or f(depvar)=-exp(depvar*c) (if c<0) instead of using the cost-function approach. This follows the approach of Bond & Lang (2019) and Kaiser & Vendrik (2023).

{p 4 4} {opt pvalue} displays p-value statistics including original p-values, minimum and maximum p-values achievable through transformations, and minimum costs needed to change statistical significance. When {cmd:pythonno} is not specified this currently only works after running {cmd:reg} and only for standard and 'robust' standard errors (and {bf:not} for clustered SEs). Also does {bf:not} work with factor variables.

{p 4 4} {cmd:critval(}{it:real}{cmd:)} sets the significance level for statistical tests. Default is 0.05. Only relevant when {cmd:pvalue} is specified. 

{p 4 4} {cmd:revpoint(}{it:real}{cmd:)} specifies the target value for coefficient reversal. Default is 0 (sign reversal).
For example, {cmd:revpoint(0.5)} checks if coefficients can be transformed to equal 0.5, and, if so, at what cost.

{marker opt_cost}{...}
{dlgtab:Cost-function options}

{p 4 4} {cmd:alpha(}{it:real}{cmd:)} controls the cost function sensitivity when using Python optimization. Higher values penalize uneven transformations more heavily.
The cost function is: cost = (var/maxvar)^(1/alpha). Default is alpha=2.

{p 4 4} {opt theil} uses the normalized Theil inequality index as the cost function instead of the alpha-based variance cost function when using Python optimization.

{marker opt_search}{...}
{dlgtab:Exponential function search options}

{p 4 4} {cmd:start(}{it:real}{cmd:)} and {cmd:end(}{it:real}{cmd:)} define the search range for the transformation parameter c in exponential transformations f(y) = exp(c*y) when using {opt pythonno}.
Wider ranges allow more extreme transformations but increase computation time.

{p 4 4} {cmd:precision(}{it:real}{cmd:)} controls the grid density for searching transformation parameters when using {opt pythonno}.
Smaller values provide more precise results but require longer computation time.

{marker opt_output}{...}
{dlgtab:Output options}

{p 4 4} {cmd:keep(}{it:string}{cmd:)} restricts the output to specific variables from the original regression.

{marker examples}{...}
{title:Examples}

{p 4 4}Basic reversal analysis after linear regression:{p_end}
{p 8 12}{inp:. sysuse auto}{p_end}
{p 8 12}{inp:. regress rep78 price mpg headroom displacement gear_ratio}{p_end}
{p 8 12}{inp:. coeff_reverser}{p_end}

{p 4 4}Include p-value bounds:{p_end}
{p 8 12}{inp:. coeff_reverser, pvalue}{p_end}

{p 4 4}Use Python optimization with Theil cost function:{p_end}
{p 8 12}{inp:. coeff_reverser, theil}{p_end}

{p 4 4}Use exponential transformation search when Python unavailable:{p_end}
{p 8 12}{inp:. coeff_reverser, pythonno}{p_end}

{p 4 4}Check reversal to specific target value (here =0.2):{p_end}
{p 8 12}{inp:. coeff_reverser, revpoint(0.2)}{p_end}

{p 4 4}Only display specific variables. Use custom exponential search range:{p_end}
{p 8 12}{inp:. coeff_reverser, pythonno keep(income education) start(-3) end(3) precision(0.05)}{p_end}

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

{p2col 5 20 24 2: Cost and c-value matrices}{p_end}
{synopt:{cmd:r(cost)}}transformation costs from Python optimization (Python mode only){p_end}
{synopt:{cmd:r(minc)}}minimum c-values for coefficient reversal ({opt pythonno} mode only){p_end}
{synopt:{cmd:r(mincp)}}minimum c-values for significance reversal ({opt pythonno} mode only){p_end}

{p2col 5 20 24 2: P-value matrices (if {opt pvalue} specified)}{p_end}
{synopt:{cmd:r(p)}}original p-values from fitted model{p_end}
{synopt:{cmd:r(minp)}}minimum p-values across transformations{p_end}
{synopt:{cmd:r(maxp)}}maximum p-values across transformations{p_end}
{synopt:{cmd:r(costp)}}transformation costs for significance reversal (Python mode only){p_end}

{p2col 5 20 24 2: Advanced matrices}{p_end}
{synopt:{cmd:r(d)}}reversal indicators for each coefficient{p_end}
{synopt:{cmd:r(hdp)}}p-values from hd transformations{p_end}
{synopt:{cmd:r(b_full)}}full coefficient matrix across transformations ({opt pythonno} mode){p_end}
{synopt:{cmd:r(r_full)}}full reversal results matrix ({opt pythonno} mode){p_end}
{synopt:{cmd:r(p_full)}}full p-value matrix ({opt pythonno} with {opt pvalue}){p_end}
{synopt:{cmd:r(y_full)}}full significance results matrix ({opt pythonno} with {opt pvalue}){p_end}

{title:References}

{p 4 4} Bond, Timothy N. and Kevin Lang. (2019). The sad truth about happiness scales. {it:Journal of Political Economy}, 127:1629–1640.

{p 4 4} Kaiser, Caspar and Anthony Lepinteur. (2025). Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data. Available at {browse "https://arxiv.org/abs/2507.16440v1":https://arxiv.org/abs/2507.16440}.

{p 4 4} Kaiser, Caspar and Maarten Vendrik. (2023). How much can we learn from happiness data? University of Oxford, Mimeo.

{title:Author}

{p 4 4}  Caspar Kaiser {p_end}
{p 4 4}  Warwick Business School, University of Warwick {p_end}
{p 4 4}  caspar.kaiser@wbs.ac.uk {p_end}

{p 4 4}  Anthony Lepinteur {p_end}
{p 4 4}  University of Luxembourg {p_end}
{p 4 4}  anthony.lepinteur@uni.lu {p_end}

{p 4 4} {cmd:coeff_reverser} is under active development. Please report bugs or feature requests to the authors.