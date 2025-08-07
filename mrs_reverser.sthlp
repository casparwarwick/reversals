{smcl}
help for {hi:mrs_reverser} version 0.1 {right: (Caspar Kaiser & Anthony Lepinteur)}
{hline}
{title:Marginal Rate of Substitution (MRS) Reversal Analysis for Ordered Response Data}

{title:Syntax}

{p 4 4 2}
{cmd:mrs_reverser}{cmd:,} {cmd:denom(}{it:varname}{cmd:)}
[{cmd:}{it:options}]

{synoptset tabbed}{...}
{synopthdr}
{synoptline}
{syntab:General {help mrs_reverser##opt_general:[+]}}
{synopt:{cmd:denom(}{it:varname}{cmd:)}}Specifies the denominator variable for ratio calculations (required){p_end}
{synopt:{opt pythonno}}Use exponential transformations instead of Python cost minimization{p_end}
{synopt:{cmd:target_ratio(}{it:real}{cmd:)}}Specifies a target ratio value to calculate transformation cost{p_end}

{syntab:Cost-function options {help mrs_reverser##opt_cost:[+]}}
{synopt:{cmd:alpha(}{it:real}{cmd:)}}Specifies the alpha parameter for the cost function (default: 2){p_end}
{synopt:{opt theil}}Use normalized Theil index as cost function (overrides {cmd:alpha} option){p_end}

{syntab:Exponential function search options (applies when specifying {cmd:pythonno}) {help mrs_reverser##opt_search:[+]}}
{synopt:{cmd:start(}{it:real}{cmd:)}}Smallest value of c over which to search (default: -2){p_end}
{synopt:{cmd:end(}{it:real}{cmd:)}}Largest value of c over which to search (default: 2){p_end}
{synopt:{cmd:precision(}{it:real}{cmd:)}}Grid precision for c values (default: 0.1){p_end}

{syntab:Output options {help mrs_reverser##opt_output:[+]}}
{synopt:{cmd:keep(}{it:string}{cmd:)}}Specifies list of variables to keep in displayed results table{p_end}

{synoptline}
{p 4 4} {cmd:mrs_reverser} requires that a regression has been run and is stored in {cmd:e()}. {p_end}
{p 4 4} Unless pythonno is specified, this command requires Python to be available, and SciPy and NumPy to be available in Stata's Python environment. {p_end}

{marker introduction}{...}
{title:Introduction}

{p 4 4}{cmd:mrs_reverser} analyzes the stability of coefficient ratios (marginal rates of substitution) to monotonic transformations of ordinal dependent variables.
This implements the methodology developed in Kaiser & Lepinteur (2025) for testing the robustness of relative effect sizes to departures from the linearity assumption in survey scales.

{p 4 4}See Kaiser & Lepinteur (2025; {browse "https://arxiv.org/abs/2507.16440v1":https://arxiv.org/abs/2507.16440v1}) for details.{p_end}

{marker description}{...}
{title:Description}

{p 4 4}{cmd:mrs_reverser} takes a linear regression and determines how coefficient ratios (marginal rates of substitution) can change through monotonic transformations of the dependent variable.

The command requires specifying a denominator variable and provides several outputs:

{p 6 8}• Original coefficient ratios (numerator/denominator) for all variables{p_end}
{p 6 8}• Minimum and maximum achievable ratios under any monotonic transformation{p_end}
{p 6 8}• Transformation cost required to achieve a specific target ratio (unless {cmd:pythonno} is specified){p_end}

Two computational approaches are available:

{p 6 8}• {bf:Python optimization}: Uses Python to find transformations minimizing a cost function as described in Kaiser & Lepinteur (2025). This finds 'least non-linear' transformations needed to achieve specific ratio targets. {p_end}
{p 6 8}• {bf:Exponential transformation search}: Should Python be unavailable, searches over exponential transformations of the form f(depvar)=exp(depvar*c) (for positive c) or f(depvar)=-exp(depvar*c) (for negative c). This follows the approach of Bond & Lang (2019) and Kaiser & Vendrik (2023). This approach is much more restrictive and only finds least non-linear reversals within the exponential class.{p_end}

{title:Options}
{marker opt_general}{...}
{dlgtab:General}

{p 4 4} {cmd:denom(}{it:varname}{cmd:)} specifies the denominator variable for all ratio calculations. This option is required.
The command will calculate ratios of all other coefficients relative to this denominator coefficient.

{p 4 4} {opt pythonno} for use when Python is not available. Searches over exponential transformations of the form f(depvar)=exp(depvar*c) (if c>0) or f(depvar)=-exp(depvar*c) (if c<0) instead of using the cost-function approach. This follows the approach of Bond & Lang (2019) and Kaiser & Vendrik (2023).

{p 4 4} {cmd:target_ratio(}{it:real}{cmd:)} specifies a target ratio value for cost calculation.
When specified, the command calculates the minimum transformation cost needed to achieve this target ratio for each numerator variable.
For example, {cmd:target_ratio(1)} calculates the cost to make each coefficient equal to the denominator coefficient. Finds smallest required "c" in exponential transformations when {cmd:pythonno} is specified. 

{marker opt_cost}{...}
{dlgtab:Cost-function options}

{p 4 4} {cmd:alpha(}{it:real}{cmd:)} controls the cost function sensitivity in Python optimization. Higher values penalize uneven transformations more heavily.
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

{p 4 4}Basic MRS analysis after linear regression:{p_end}
{p 8 12}{inp:. sysuse auto}{p_end}
{p 8 12}{inp:. regress rep78 price mpg headroom displacement gear_ratio}{p_end}
{p 8 12}{inp:. mrs_reverser, denom(price)}{p_end}

{p 4 4}Calculate cost to achieve specific target ratio:{p_end}
{p 8 12}{inp:. mrs_reverser, denom(mpg) target_ratio(0.5)}{p_end}

{p 4 4}Use Theil cost function:{p_end}
{p 8 12}{inp:. mrs_reverser, denom(headroom) theil}{p_end}

{p 4 4}Use exponential transformation search when Python unavailable:{p_end}
{p 8 12}{inp:. mrs_reverser, denom(displacement) pythonno}{p_end}

{p 4 4}Focus on specific variables with custom alpha parameter:{p_end}
{p 8 12}{inp:. mrs_reverser, denom(displacement) keep(price mpg) alpha(3)}{p_end}

{p 4 4}Calculate cost for ratio equal to 1 (numerator equals denominator):{p_end}
{p 8 12}{inp:. mrs_reverser, denom(gear_ratio) target_ratio(1) alpha(1.5)}{p_end}

{p 4 4}Exponential search with custom range and precision:{p_end}
{p 8 12}{inp:. mrs_reverser, denom(price) pythonno start(-3) end(3) precision(0.05)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mrs_reverser} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Main results}{p_end}
{synopt:{cmd:r(result)}}complete displayed table with all computed statistics{p_end}

{p2col 5 20 24 2: Coefficient ratio matrices}{p_end}
{synopt:{cmd:r(ratio)}}original coefficient ratios (numerator/denominator){p_end}
{synopt:{cmd:r(minratio)}}lower bounds for coefficient ratios{p_end}
{synopt:{cmd:r(maxratio)}}upper bounds for coefficient ratios{p_end}

{p2col 5 20 24 2: Cost matrices (if {opt target_ratio} specified)}{p_end}
{synopt:{cmd:r(cost)}}transformation costs to achieve target ratio (Python mode only){p_end}
{synopt:{cmd:r(minc)}}minimum c-values for achieving target ratio ({opt pythonno} mode only){p_end}

{marker technical}{...}
{title:Technical notes}

{p 4 4} Before performing the analysis, the command checks whether the denominator coefficient can be sign-reversed through transformations.
If the denominator is reversible, ratio bounds become infinite. A warning is displayed in that case as results may become unreliable in that case.

{title:References}

{p 4 4} Bond, Timothy N. and Kevin Lang. (2019). The sad truth about happiness scales. {it:Journal of Political Economy}, 127:1629–1640.

{p 4 4} Kaiser, Caspar and Anthony Lepinteur. (2025). Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data. Available at {browse "https://arxiv.org/abs/2507.16440v1":https://arxiv.org/abs/2507.16440v1}.

{p 4 4} Kaiser, Caspar and Maarten Vendrik. (2023). How much can we learn from happiness data? University of Oxford, Mimeo.

{title:Author}

{p 4 4}  Caspar Kaiser {p_end}
{p 4 4}  Warwick Business School, University of Warwick {p_end}
{p 4 4}  caspar.kaiser@wbs.ac.uk {p_end}

{p 4 4}  Anthony Lepinteur {p_end}
{p 4 4}  University of Luxembourg {p_end}
{p 4 4}  anthony.lepinteur@uni.lu {p_end}

{p 4 4} {cmd:mrs_reverser} is under active development. Please report bugs or feature requests to the authors.