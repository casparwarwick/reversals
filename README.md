# REVERSALS: Robustness Analyses of Scale Transformations in Linear Regressions

This Stata package implements the methodology from Kaiser & Lepinteur (2025; https://arxiv.org/abs/2507.16440).

## Installation

### From GitHub
```stata
net install reversals, from(https://raw.githubusercontent.com/casparwarwick/reversals/main/)
```

You can uninstall the package from Stata by running:

```stata
ado uninstall reversals
```

## Requirements

### Python

Python (including NumPy and SciPy) needs to be available in Stata for full functionality. 

You can check the availability of Python in Stata by typing:

```stata
python query
```
If Python and/or the NumPy and SciPy libraries are unavailable, both `coeff_reverser` and `mrs_reverser` will revert to `pythonno`.

If Numpy and/or SciPy are not installed to your Python environment, you should be able to install them by typing `pip install numpy` and `pip install scipy` in the terminal (on macOS) or the command prompt (on Windows).

If you get stuck do feel free to contact me (caspar.kaiser@wbs.ac.uk). 

### Stata

Tested on Stata 17.0. 

Requires:
- `estout` by Ben Jann. Can be installed via `ssc install estout`.
- `matselrc` (from "Yet more matrix commands") by Nick Cox. Can be installed via `net install dm79, from(http://www.stata.com/stb/stb56)`.

## Commands

### `coeff_reverser`
Tests coefficient robustness to scale transformations using one of two approaches:

- **Python optimization** (default): Uses Python to find minimal-cost transformations via numerical optimization.
- **Exponential search** (`pythonno` option): Searches over exponential transformations of the form f(depvar)=exp(c*depvar).

### `mrs_reverser` 
Tests coefficient ratio stability to scale transformations:

- Requires specifying one denominator variable via the `denom(varlist)` option.
- Returns original ratio, min/max bounds, and cost to achieve target ratio.

## Quick Start Example

```stata
sysuse auto, clear

* Run regression
reg price mpg weight length turn gear_ratio

* Test coefficient robustness
coeff_reverser

* Test ratio robustness  
mrs_reverser, denom(weight)
```

## Files

## Core Commands
- `coeff_reverser.ado`: Main coefficient reversal command.
- `mrs_reverser.ado`: Main coefficient ratio (MRS) analysis command.

## Help Files
- `coeff_reverser.sthlp`
- `mrs_reverser.sthlp`

## Python Optimization Scripts  
- `sign_reversal_cost_minimizer.py`: Coefficient reversal analysis using SciPy.
- `p_value_cost_minimizer.py`:  P-value analysis using SciPy.
- `mrs_reverser_python.py`: Coefficient ratio analysis using SciPy.

## Citation

If you use this package, please cite:

Kaiser, C. & Lepinteur, A. (2025). "Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data." https://arxiv.org/abs/2507.16440

## Author

**Caspar Kaiser**  
Warwick Business School  
Email: caspar.kaiser@wbs.ac.uk

## Support

For questions, bug reports, or feature requests, please:
1. Check the help files: `help coeff_reverser` and `help mrs_reverser`
2. Email: caspar.kaiser@wbs.ac.uk
