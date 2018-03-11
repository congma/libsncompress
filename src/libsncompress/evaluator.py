# vim: spell spelllang=en
import numpy
import scipy.linalg
import scipy.optimize
from . import simplecache


def _initial_guess(base):
    # Guess value for the three standardization parameters.
    nuis_guess = numpy.array((0.13, -3, 0.0))
    iids = numpy.fromiter(base.binidcontent[0], dtype=int)
    fids = numpy.fromiter(base.binidcontent[base.bins.nbins - 1], dtype=int)
    xi, xf = base.logredshifts[iids], base.logredshifts[fids]
    yi, yf = base.peakbmags[iids], base.peakbmags[fids]
    xim = xi.mean()
    yim = yi.mean()
    if base.bins.nbins == 1:
        slope = 0.0
    else:
        slope = (yf.mean() - yim) / (xf.mean() - xim)
    dist_guess = numpy.array([yim + (lz - xim) * slope for lz in
                              base.bins.itercontrolpoints()]) + 19.05
    return numpy.concatenate((nuis_guess, dist_guess))


def _scale_guess(cplist):
    nuis_guess = numpy.array((0.1, 1.3, 0.4))
    width = numpy.log10(1.3) + 2.0
    x = numpy.array([(v + 2.0) / width for v in cplist])
    l2 = numpy.log(2.0)
    l4 = l2 * 2.0
    k1 = -l4 / 0.618
    k2 = l4 / 0.382

    @numpy.vectorize
    def intp(t):
        if t < 0.618:
            y = l2 + k1 * t
        else:
            y = -l2 + k2 * (t - 0.618)
        return y

    dist_guess = numpy.exp(intp(x))
    return numpy.concatenate((nuis_guess, dist_guess))


def lmultiply_inv(a_factor, b):
    """Compute the multiplication A^-1 B given the Cholesky factor of A."""
    return scipy.linalg.cho_solve(a_factor, b, check_finite=False)


def quadratic_inv(matrix_factor, lvector, rvector):
    """Compute the quadratic form < x, A^-1 y > given lvector
    as x and rvector as y.  The input matrix_factor is the Cholesky factor
    of A.
    """
    rcomplement = lmultiply_inv(matrix_factor, rvector)
    return numpy.dot(lvector, rcomplement)


def logdet_cholesky(matrix_factor):
    """Compute the logarithm of determinant of matrix given by its Cholesky
    factor.  Only for positive-determined matrices.
    """
    return 2.0 * numpy.sum(numpy.log(numpy.diag(matrix_factor[0])))


def traceofdot(a, b):
    """Equivalent to numpy.trace(numpy.dot(a, b)), but use einsum so possibly
    faster than the naive implementation.
    """
    return numpy.einsum("ij, ji", a, b)


def hesstocov(matrix):
    """Convert Hessian matrix to covariance matrix estimate."""
    hessfac = scipy.linalg.cho_factor(matrix, check_finite=False)
    covh = lmultiply_inv(hessfac, numpy.eye(matrix.shape[0]))
    return covh + covh.T


class CovEvaluator(simplecache.ArrayMethodCacheMixin, object):
    """Likelihood evaluation for the binned supernovae.
    We follow the convention that the residual vector is data minus
    "model", not the other way around.

    Attributes:
        dimension:  number of parameters
        base:  BinnedSN instance
        withlogdet:  boolean, whether the correct conditional probability is
                     used
        res:  evaluator output; initialized to None
    """
    def __init__(self, basesn, withlogdet=True, *args, **kwargs):
        """Set-up internal storage of covariance evaluator.
        Arguments:
            basesn:  BinnedSN instance
            withlogdet:  boolean, whether the correct conditional probability
                         is used
            *args, **kwargs will be passed to the superclass initializer
        """
        super(CovEvaluator, self).__init__(*args, **kwargs)
        self.res = None
        self.dimension = 3 + basesn.bins.ncontrolpoints  # n of parameters
        self.base = basesn
        self.withlogdet = withlogdet
        # Gradient of residual vector.
        self._grad_vec = (self._mk_grad_datavec_const() -
                          self._mk_grad_fitvec_const())
        # 2nd derivative of the full, data covariance matrix
        # with respect to alpha and beta
        self._hesscov = self._mk_hess_cov_const()
        return None

    def _grad_cov(self, param):
        """Calculate the partial derivatives of the covariance matrix with
        respect to alpha and beta.
        Returns a 2-array containing the partial derivative matrices.
        """
        res = numpy.empty((2,
                           self.base.datadimension, self.base.datadimension))
        coeffs = numpy.array((1.0, param[0], param[1]))
        for k in (1, 2):
            tmp = numpy.einsum("i, ijl",
                               coeffs, self.base.bblock[:, k, :, :])
            tmp += numpy.einsum("i, ijl",
                                coeffs, self.base.bblock[k, :, :, :])
            res[k - 1] = tmp
        return res

    def _fill_deriv_bins(self, vec, binid, left=True):
        """Fill the array vec with derivatives of binned mu's for the given
        bin id and a boolean flag indicating whether this bin is considered
        to the left or right of the control point.

        vec is modified in-place.  Returns None.
        """
        thisset = self.base.binidcontent[binid]
        lower, upper = self.base.bins.getbin(binid)
        rx = upper - lower
        for dataid in thisset:
            lz = self.base.logredshifts[dataid]
            sl = (lz - lower) / rx
            vec[dataid] = sl if left else (1.0 - sl)
        return None

    def _mk_hess_cov_const(self):
        """Similar to _grad_cov(), but for the 2nd-order derivatives, or
        Hessian.
        """
        res = [[], []]
        indices = (1, 2)
        for i in indices:
            for j in indices:
                res[i - 1].append(self.base.bblock[i][j] +
                                  self.base.bblock[j][i])
        return numpy.ascontiguousarray(res)

    def _mk_grad_datavec_const(self):
        """Calculate the derivative of datavec(param).
        Each row of the returned array is the derivative with respect
        to the corresponding parameter.
        """
        res = numpy.zeros((self.dimension, self.base.datadimension))
        # alpha, beta
        res[0, :] = self.base.shapes
        res[1, :] = self.base.colors
        # delta
        res[2, :] = numpy.where(self.base.needhostcorr, -1.0, 0.0)
        # the mu_b's needs no further computing.
        return res

    def _mk_grad_fitvec_const(self):
        """Like _mk_grad_datavec_const, but for the fit vector."""
        res = numpy.zeros((self.dimension, self.base.datadimension))
        # alpha, beta and delta needs no further computing.
        # the mu_b's
        row = 3
        for (endpoint,
             leftbin, rightbin) in self.base.bins.itercontrolpoints_classify():
            if leftbin is not None:
                self._fill_deriv_bins(res[row, :], leftbin, left=True)
            if rightbin is not None:
                self._fill_deriv_bins(res[row, :], rightbin, left=False)
            row += 1
        return res

    @simplecache.memoized()
    def _comp_interm_for_chisq(self, param):
        # Data (distance modulus) cov at param
        cov = self.cov(param)
        # Cholesky factor of cov, to be used for computing cov^{-1} x
        covfac = scipy.linalg.cho_factor(cov, check_finite=False)
        # Data residual vector at param
        vec = self.datavec(param) - self.fitvec(param)
        return cov, covfac, vec

    @simplecache.memoized()
    def _comp_interm_for_grad_hess(self, param):
        cov, covfac, vec = self._comp_interm_for_chisq(param)
        # rightfactor = cov^{-1} vec
        rightfactor = lmultiply_inv(covfac, vec)
        # deriv_cov = jacobian of covariance
        deriv_cov = self._grad_cov(param)
        ddr = numpy.dot(deriv_cov, rightfactor)
        return rightfactor, deriv_cov, ddr

    @simplecache.memoized()
    def _comp_interm_for_hess_corr(self, param):
        cov, covfac, vec = self._comp_interm_for_chisq(param)
        rightfactor, deriv_cov, ddr = self._comp_interm_for_grad_hess(param)
        return [lmultiply_inv(covfac, d) for d in deriv_cov]

    def cov(self, param):
        """Evaluate the covariance of data at given parameter value."""
        alpha = param[0]
        beta = param[1]
        parray = numpy.array((1.0, alpha, beta))
        datacov = numpy.einsum("ijkl, i, j", self.base.bblock, parray, parray)
        # symmetrize for stability?
        datacov = (datacov + datacov.T) / 2.0
        return datacov

    def datavec(self, param):
        """Evaluate the data vector (i.e. the distance moduli) at the parameter
        value.
        """
        alpha = param[0]
        beta = param[1]
        dmb = param[2]
        # NOTE: M_b^1 fixed to -19.05
        vec = (self.base.peakbmags +
               alpha * self.base.shapes +
               beta * self.base.colors + 19.05)
        # Host correction, mass-dependent term
        vec[self.base.needhostcorr] -= dmb
        return vec

    def fitvec(self, param):
        """Evaluate the piecewise-linear fit to distance moduli.
        Returns a vector of reconstructed distance moduli evaluated at the
        given parameter values.
        """
        mus = param[3:]
        rec = numpy.empty(self.base.datadimension)
        for i, lz in enumerate(self.base.logredshifts):
            nbin, (lower, upper) = self.base.binnings[i]
            slope = (lz - lower) / (upper - lower)
            rec[i] = ((1.0 - slope) * mus[nbin] + slope * mus[nbin + 1])
        return rec

    def chisq(self, param):
        """Evaluate the "effective chi-squared" statistic at the parameter
        value.  It is actually -2 * log-likelihood without the constant term
        len(data) * log (2 * pi).  If self.withlogdet is False, it evaluates to
        the literal "chi-squared".
        """
        cov, covfac, vec = self._comp_interm_for_chisq(param)
        res = quadratic_inv(covfac, vec, vec)
        if self.withlogdet:
            res += logdet_cholesky(covfac)
        return res

    def chisqgrad(self, param):
        """Evaluate the derivative of effective chisq."""
        cov, covfac, vec = self._comp_interm_for_chisq(param)
        rightfactor, deriv_cov, ddr = self._comp_interm_for_grad_hess(param)
        if self.withlogdet:
            cc = self._comp_interm_for_hess_corr(param)
        res = 2.0 * numpy.dot(self._grad_vec, rightfactor)
        # NOTE: Still a small-scale explicit for-loop
        for i in (0, 1):
            res[i] -= numpy.dot(ddr[i], rightfactor)
            if self.withlogdet:
                res[i] += numpy.trace(cc[i])
        return res

    def chisqhess(self, param):
        """Evaluate the Hessian of effective chisq."""
        cov, covfac, vec = self._comp_interm_for_chisq(param)
        rightfactor, deriv_cov, ddr = self._comp_interm_for_grad_hess(param)
        if self.withlogdet:
            cc = self._comp_interm_for_hess_corr(param)
        # The following line does the magic
        res = 2.0 * quadratic_inv(covfac, self._grad_vec, self._grad_vec.T)
        # Partially "unrolled" loops
        for i, j in ((0, 0), (0, 1), (1, 1)):
            pi = ddr[i]
            pj = ddr[j]
            # This is manifestly symmetric in i, j
            cij = self._hesscov[i][j]
            # Symmetric, too.
            res[i][j] += 2.0 * quadratic_inv(covfac, pj, pi)
            # Sum of the following two terms are symmetric.
            res[i][j] -= 2.0 * numpy.dot(lmultiply_inv(covfac,
                                                       self._grad_vec[i]),
                                         pj)
            res[i][j] -= 2.0 * numpy.dot(lmultiply_inv(covfac,
                                                       self._grad_vec[j]),
                                         pi)
            # Manifestly symmetric.
            res[i][j] -= numpy.dot(rightfactor, numpy.dot(cij,
                                                          rightfactor))
            # NOTE: This code block is the greatest hotspot of all.
            # Any optimization would be much much appreciated.
            if self.withlogdet:
                # Trace of dot is symmetric with the two matrices being dotted.
                res[i][j] -= traceofdot(cc[j], cc[i])
                # cij is manifestly symmetric
                res[i][j] += numpy.trace(lmultiply_inv(covfac, cij))
        res[1][0] = res[0][1]
        return res

    def minimize(self, x0=None, jac=None, hess=None,
                 xscalings=None, fscaling=1.0, *otherargs, **kwargs):
        """Wrapper around scipy.optimize.minimize, with scaling.

        Output is scaled back using scaling inputs.  If not given, use
        guesstimates. When properly chosen, scaling as a form of conditioning
        may improve convergence speed.

        ``x0`` is the initial guess, and if not given, or is None, defaults to
        a guesstimate.

        ``jac`` and ``hess`` are functions that evaluates the Jacobian and
        Hessian respectively.  If not given or are None, use the default ones
        supplied by self.

        Other positional arguments ``otherargs`` are passed to
        scipy.optimize.minimize.  Other keyword arguments ``kwargs`` are also
        passed through, but if the ``method`` keyword is absent, use
        ``trust-ncg`` by default.

        Returns the scipy.optimize.OptimizeResult object and set it to the
        ``res`` attribute of self instance.  The OptimizeResult object has an
        additional attribute, ``compressed_cov``, obtained by casting the
        Hessian at optimal point as an estimate of the covariance matrix of
        compressed results.
        """
        newkw = dict(kwargs, method=kwargs.get("method", "trust-ncg"))
        if x0 is None:
            x0 = _initial_guess(self.base)
        if xscalings is None:
            local_scalings = _scale_guess(self.base._cpcache)
        else:
            local_scalings = xscalings
        assert local_scalings.shape == x0.shape
        assert not numpy.allclose(local_scalings, 0)
        assert fscaling > 0
        if jac is None:
            thisjac = self.chisqgrad
        else:
            thisjac = jac
        if hess is None:
            thishess = self.chisqhess
        else:
            thishess = hess
        init = x0 / local_scalings
        outer = numpy.outer(local_scalings, local_scalings)
        # Optimize by pre-computing constants, at the cost of higher memory
        # footprint.
        fs_local = fscaling * local_scalings
        fs_outer = fscaling * outer

        def funwrap(scaled_x):
            return fscaling * self.chisq(scaled_x * local_scalings)

        def gradwrap(scaled_x):
            return fs_local * thisjac(scaled_x * local_scalings)

        def hesswrap(scaled_x):
            return fs_outer * thishess(scaled_x * local_scalings)

        res = scipy.optimize.minimize(funwrap, init,
                                      jac=gradwrap, hess=hesswrap,
                                      *otherargs, **newkw)
        res.x *= local_scalings
        res.fun /= fscaling
        try:
            res.jac /= fs_local
        except AttributeError:
            pass
        try:
            res.hess /= fs_outer
        except AttributeError:
            pass
        # Do not retrieve hess from res; as noted above, not all minimizers
        # support it.  chisqhess method should be fast enough due to caching.
        res.compressed_cov = hesstocov(self.chisqhess(res.x))
        self.res = res
        return res

    @property
    def compressed_cov(self):
        """Convenient property to access the covariance estimate from
        self.minimize() run, if it is available.
        """
        if self.res is None or not self.res.success:
            raise ValueError("Minimization result unusable.")
        return self.res.compressed_cov
