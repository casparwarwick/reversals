# REVERSALS: Robustness Analyses of Scale Transformations in Linear Regressions

This Stata package implements the methodology from Kaiser & Lepinteur (2025; https://arxiv.org/abs/2507.16440v1).

## Installation

### From GitHub
```stata
net install reversals, from(https://raw.githubusercontent.com/casparwarwick/reversals/main/)
```

Note that python (including the SciPy and NumPy packages) needs to be available in Stata (unless you specifiy `pythonno`).

You can check the availability of Python in Stata by typing:

```stata
python query
```
Do feel free to contact me (caspar.kaiser@wbs.ac.uk) in case you get stuck. 

You can install the package from Stata by running:

```stata
ado uninstall reversals
```

## Commands

### `coeff_reverser`
Tests coefficient robustness to scale transformations using one of two approaches:

- **Python optimization** (default): Uses Python to find minimal-cost transformations via numerical optimization.
- **Exponential search** (`pythonno` option): Searches over exponential transformations f(y)=exp(c*y).

### `mrs_reverser` 
Tests coefficient ratio stability to scale transformations:

- Requires specifying one denominator variable via `denom(varlist)` option
- Returns original ratio, min/max bounds, and cost to achieve target ratio

## Requirements

- **Stata**: Tested on Version 17.0 (requires Python integration)
- **Python** with packages:
  - `numpy`
  - `scipy` 
  - `sfi` (Stata Function Interface)

## Quick Start Example

```stata
sysuse auto, clear

* Run regression
reg price mpg weight

* Test coefficient robustness
coeff_reverser

* Test ratio robustness  
mrs_reverser, denom(weight)
```

## Files

## Core Commands
- `coeff_reverser.ado`: Main coefficient reversal command
- `mrs_reverser.ado`: Main coefficient ratio (MRS) analysis command  

## Help Files
- `coeff_reverser.sthlp`
- `mrs_reverser.sthlp`

## Python Optimization Scripts  
- `sign_reversal_cost_minimizer.py`: Coefficient reversal analysis using SciPy
- `p_value_cost_minimizer.py`:  p-value analysis using SciPy
- `mrs_reverser_python.py`: Coefficient ratio analysis using SciPy

## Citation

If you use this package, please cite:

Kaiser, C. & Lepinteur, A. (2025). "Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data." https://arxiv.org/abs/2507.16440v1

## Author

**Caspar Kaiser**  
Warwick Business School  
Email: caspar.kaiser@wbs.ac.uk

## Support

For questions, bug reports, or feature requests, please:
1. Check the help files: `help coeff_reverser` and `help mrs_reverser`
2. Email: caspar.kaiser@wbs.ac.uk
