## CalcGamma

by Mark Geurts <mark.w.geurts@gmail.com>
<br>Copyright &copy; 2015, University of Wisconsin Board of Regents

## Description

`CalcGamma()` computes a 1-D, 2-D, or 3-D local or global Gamma index between two datasets (reference and target) given a defined coordinate space using MATLAB.  The Gamma analysis is performed based on the formalism presented by D. A. Low et. al., [A technique for the quantitative evaluation of dose distributions.](http://www.ncbi.nlm.nih.gov/pubmed/9608475), Med Phys. 1998 May; 25(5): 656-61.

## Installation

To install this function, copy `CalcGamma.m` from this repository into your MATLAB path. If installing as a submodule into another git repository, execute `git submodule add https://github.com/mwgeurts/gamma`.  

## Usage and Documentation

```matlab
gamma = CalcGamma(reference, target, percent, dta);
gamma = CalcGamma(..., 'OptionName', OptionValue);
```

See the [wiki](../../wiki) for information on the format of the input and arguments as well options and examples.

## License

Released under the GNU GPL v3.0 License.  See the [LICENSE](LICENSE) file for further details.
