# REVERSALS: Robust Analysis of Scale Transformations in Ordinal Data

This Stata package implements the methodology from Kaiser & Lepinteur (2025) "Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data" for analyzing the robustness of empirical results to departures from linearity assumptions in ordinal scales.

## Installation

### From GitHub
```stata
net install reversals, from(https://raw.githubusercontent.com/[username]/reversals/main/)
```

### Manual Installation
1. Download all files from this repository
2. Place them in your Stata ado-path
3. Ensure Python integration is available in Stata

## Commands

### `coeff_reverser`
Tests coefficient robustness to scale transformations using two computational approaches:

- **Python optimization** (default): Uses SciPy to find minimal-cost transformations via global optimization
- **Exponential search** (`pythonno` option): Searches over exponential transformations f(y)=exp(c*y) following Bond & Lang (2019)

**Key Features:**
- Returns coefficient bounds, transformation costs, and minimum c-values for significance changes  
- Advanced p-value bounds analysis with cost minimization
- Supports multiple regression types and cost functions
- Comprehensive stored results in `r()` matrices

**Basic Syntax:**
```stata
coeff_reverser varlist [if] [in], [options]
```

### `mrs_reverser` 
Tests coefficient ratio stability to scale transformations:

- Requires specifying one denominator variable via `denom(varlist)` option
- Checks if denominator is reversible (error if yes)  
- Returns original ratio, min/max bounds, and cost to achieve target ratio
- Supports `alpha()` and `theil` cost function options

**Basic Syntax:**
```stata  
mrs_reverser varlist [if] [in], denom(varlist) [options]
```

## Requirements

- **Stata**: Version 16+ (requires Python integration)
- **Python**: Version 3.6+ with packages:
  - `numpy`
  - `scipy` 
  - `sfi` (Stata Function Interface)

## Quick Start Example

```stata
sysuse auto, clear

* Test coefficient robustness
coeff_reverser price mpg weight, pvalue

* Test ratio robustness  
mrs_reverser price mpg, denom(weight)
```

## Methodology

The framework applies to any bounded ordered scale (wellbeing, risk preferences, trust, etc.) and addresses the fundamental question of whether empirical results depend on implicit linearity assumptions in ordinal data.

**Key Concepts:**
- **Cost function** C ∈ [0,1]: Measures deviation from linearity (0=perfectly linear, 1=maximally non-linear)
- **Scale transformations**: Monotonic functions that preserve ordinality while changing implied cardinality
- **Coefficient bounds**: Min/max achievable coefficients under cost constraints
- **Ratio stability**: Robustness of coefficient ratios (e.g., marginal rates of substitution)

## Files

### Core Commands
- `coeff_reverser.do`: Main coefficient reversal analysis command
- `mrs_reverser.do`: Main coefficient ratio (MRS) reversal analysis command  
- `coeff_reverser.sthlp`: Comprehensive help file with examples
- `mrs_reverser.sthlp`: Comprehensive help file for MRS analysis

### Python Optimization Scripts  
- `sign_reversal_cost_minimizer.py`: Coefficient reversal optimization using SciPy
- `p_value_cost_minimizer.py`: Advanced p-value bounds analysis with dual optimization
- `mrs_reverser_python.py`: Coefficient ratio optimization using SciPy

### Package Files
- `reversals.pkg`: Stata package description file
- `stata.toc`: Table of contents for Stata installation
- `LICENSE`: MIT License

## Citation

If you use this package, please cite:

Kaiser, C. & Lepinteur, A. (2025). "Measuring the unmeasurable? Systematic evidence on scale transformations in subjective survey data." *Journal of Applied Econometrics* (forthcoming).

## Author

**Caspar Kaiser**  
Warwick Business School  
Email: caspar.kaiser@wbs.ac.uk

## License

MIT License - see LICENSE file for details.

## Support

For questions, bug reports, or feature requests, please:
1. Check the help files: `help coeff_reverser` and `help mrs_reverser`
2. Email: caspar.kaiser@wbs.ac.uk
3. Open an issue on GitHub

## Version History

- **v1.0.0** (August 2025): Initial release with complete `coeff_reverser` and `mrs_reverser` implementation
- Full Python-Stata integration via SFI
- Comprehensive p-value bounds analysis  
- Advanced optimization algorithms with constraint handling