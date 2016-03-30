## Summary ##

`libsncompress` -- a Python utility for compressing supernova cosmological 
data.

## Description ##

The Python package `libsncompress` implements the linear compression method
described in the paper [Ma et al. (2016)][m16].  It is designed for use with
the [JLA][jla] dataset, but can be easily extended for other similar datasets.

It also includes a Python executable script, 
[`jlacompress`](scripts/jlacompress), that serves as an example command-line 
user interface.

## JLA Compression Script ##

### Synopsis ###

```
jlacompress [-h] [-d DIR] [-t FILE] [-p PREFIX] [-c z [z ...]] [-n] [-v]
```

The script requires JLA data files to run.  See the section [Data 
Files](#data-files) for details.

### Options ###

*  `-h`, `--help`: show help message and exit 
*  `-d DIR`, `--covdir DIR`:  path to directory of <abbr
   title="Flexible Image Transport System">FITS</abbr> covariance files 
   (default: `./covmat`)
*  `-t FILE`, `--table FILE`: path to JLA data table file (default: 
   `./jla_lcparams.txt`)
*  `-p PREFIX`, `--output-prefix PREFIX`: prefix to output file names
*  `-c z [z ...]`, `--controls z [z ...]`: locations of control points (as 
   redshift).  At least two control points are required.  If unspecified, use 
   the default control points in the [JLA paper][jla].
*  `-n`, `--no-logdet`: don't use the correct conditional probability 
   (default: off; *warning:* use at your own risk)
*  `-v`, `--verbose`: turn on verbose output (default: off)

### Output ###

The script writes three output files when it solves the optimization problem 
successfully:

*  Mean (approximate, actually posterior-optimizing) compression parameters in 
   the order of (`alpha`, `beta`, `delta_M`, `mu1`, `mu2`, ...), where `muN` 
   is the value of distance modulus at the `N`-th control point.
*  Covariance matrix (symmetric but with all elements filled) of the 
   parameters.
*  List of redshift of control points.

The file names are `mean.txt`, `cov.txt` and `redshift.txt` respectively, but 
they can be prefixed by arbitrary strings specified by the user with the 
`-p`/`--prefix` option.

The first three lines in `mean.txt` are post-compression estimates of best-fit 
parameters for standardization parameters, in the order of `alpha`, `beta` and 
`delta`.  The rest of the lines are for the compressed distance modulus at 
each control point, in the order of increasing redshift.  In `cov.txt`, the 
rows and columns are in the same order.

If the prefix contains slash (`/`) characters, it will be understood as 
directory separators.  The files will fail to write to their paths if the 
nested directory paths do not exist.

The output files will have suffix `-no-logdet` appended to the path but before 
the `.txt` extension, if `-n`/`--no-logdet` is specified.

When verbose output is enabled by `-v`/`--verbose`, additional text will be 
written to the standard error.

### Exit Status ###

The script exits with `0` for success.  Any other value indicates error.

### Example ###

Assuming the data files are in their default locations, the following command 
reproduces the default compression results in the [JLA paper][jla].

```bash
jlacompress -n
```

## Data Files ##

The JLA data files are *required* for using the package.  However, we cannot 
distribute them with the source package.  Please read the [JLA readme][jlarm]
page for details about the data files.

The following *two* files must be downloaded:

1.  The file [`jla_likelihood_v6.tgz`][jlatarball], compressed archive 
    containing the file `data/jla_lcparams.txt`.  This file contains the 
    supernova sample catalogue.  The other files in this archive are not 
    necessary.
2.  The FITS files containing the components of data covariance, in the 
    compressed archive [`covmat_v6.tgz`][jlafits].  The non-FITS files in this 
    archive are not necessary.

## Hacking ##

To use the package directly in your own Python program, simply

```python
import libsncompress
```

This will import three classes from its sub-modules into the `libsncompress` 
namespace:

*  `BinnedSN`:  data-file loader and pre-processor
*  `BinCollection`:  redshift binning and sanitizer; not very useful on its 
   own
*  `CovEvaluator`:  the actual compressor

The first thing you need to do is to specify a list (or `numpy` array) of 
control points, by their *base-10 logarithm* values.  Currently, due to 
development legacy, the "binning" class and methods are not particularly 
efficient.  This is usually not a problem because it will be used only once.

This list or array of control points must be encapsulate in *another* 
container (list, array, or tuple, etc.) before passing to the initializer of 
`libsncompress.BinnedSN` class.  The instance can be initialized by

```python
binned_sn = libsncompress.BinnedSN(basedirpath,
                                   tablepath,
                                   logbins=control_points)
```

Here `basedirpath` is the path to the directory containing the FITS covariance 
data files, `tablepath` the path to the text file containing the JLA dataset 
table, and `logbins` is the nested list of control points just obtained.

After this, we can initialize the evaluator `libsncompress.CovEvaluator` 
class, which implements the evaluation of probability log-density functions 
and their first 2 derivatives, like this:

```python
ev = libsncompress.CovEvaluator(binned_sn, withlogdet=True)
```

The optional argument `withlogdet` controls whether the full effect of 
parameter-dependent covariance matrix is taken into account.  It is so named 
due to the ubiquitous presence of "ln det Cov" term.  It defaults to `True` 
but can be set to `False`, which will evaluate the functions as if the 
customary chi-squared method were used.

The `CovEvaluator` instance, `ev`, provides a method `minimize`, which is a 
wrapper of `scipy.optimize.minimize`.  Additional positional and keyword 
arguments are passed over to that function.  The recommended optimization 
algorithm is `trust-ncg` which fully utilizes the Hessian matrix.  This can be 
enabled by passing `method="trust-ncg"` as an optional keyword parameter.

The return value of `CovEvaluator.minimize` method is simply that of the 
underlying `scipy` function, but with results suitably scaled.

The Hessian of log-PDF function can be obtained, then, at the minimizing point 
in the parameter space.  This can be used for constructing the approximate 
covariance of compression parameters.

Please notice that this implementation here is not a general, abstract 
implementation of the linear compression method detailed in [our paper][m16].
It specifically implements the sawtooth-basis compression, which is compatible 
with the original [JLA one][jla].  The implementation details, as well as the 
exposed API, are likely to see significant revisions in the future.

## Installation ##

The package source directory [`libsncompress`](libsncompress/) can be used 
directly without installation.  The package and script can also be installed 
using the standard `distutils` setup script:

```bash
python setup.py install
```

## Requirements ##

*  [`numpy`][numpy] (`>= 1.6.0`), for array data structure and basic 
   operations;
*  [`scipy`][scipy] (`>= 0.11.0`), for linear algebra and numerical 
   optimization;
*  [`pyfits`][pyfits] (unknown version), for loading FITS files;
*  [`cachetools`][ct] (unknown version), for caching partial evaluation 
   results, which is essential for compression speed.

## Performance Notes ##

Performance is mostly determined by the following two conditions:

1.  Underlying <abbr title="Basic Linear Algebra Subprograms">BLAS</abbr>/LAPACK
    libraries used by `numpy`/`scipy`, especially the "linear solver by 
    Cholesky decomposition", `(D)POTRS` function of LAPACK.  For [NetLib 
    LAPACK][netliblapack], this in turn is largely determined by the speed of 
    the level-3 BLAS triangular solver, `(D)TRSM`.  The NetLib reference 
    implementation is rather naive, and an optimized implementation of BLAS is 
    likely to boost the performance.
2.  Choice of initial value and scaling for numerical optimization.  If 
    they are suitably chosen, the number of iterations required to achieve 
    convergence is reduced.

The script [`jlacompress`](scripts/jlacompress) attempts to automatically 
create acceptable initial value and scaling that is optimized for the 
*default* compression used in the [JLA paper][jla].  The automatic initial 
value and scaling are not optimized for any other usage cases.

## Issue Tracker ##

Please report problems via the [issue tracker][issues].


[m16]: http://arxiv.org/abs/1603.08519 "The M16 paper (preprint)"
[jla]: http://arxiv.org/abs/1401.4064 "JLA reference paper"
[jlarm]: http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html "JLA project"
[jlatarball]: http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz
[jlafits]: http://supernovae.in2p3.fr/sdss_snls_jla/covmat_v6.tgz
[numpy]: http://www.numpy.org/ "NumPy homepage"
[scipy]: https://www.scipy.org/ "SciPy homepage"
[pyfits]: https://pythonhosted.org/pyfits/ "PyFITS"
[ct]: https://pythonhosted.org/cachetools/ "cachetools"
[issues]: https://gitlab.com/congma/libsncompress/issues "Issue tracker"
[netliblapack]: http://www.netlib.org/lapack/ "NetLib LAPACK"

<!--
vim: ft=markdown tw=78 fo+=tqwn spell spelllang=en et ts=4
-->
