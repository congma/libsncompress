# vim: spell spelllang=en


import os.path
import collections
import six.moves as sm
import numpy
from astropy.io import fits
from . import binning


# Table 9 of B14
_SIGMA_COH_TABLE = {1: 0.08,     # SNLS
                    2: 0.108,    # SDSS-II
                    3: 0.134,    # low-z, CfA
                    4: 0.1,      # HST
                    }
_COEFF_SIGMAZ = 1.0864878442920631136702801549045434e-3


def _preeval_diag_in_cov(redshifts, setnum):
    """Pre-compute the constant diagonal due to peculiar redshift, lensing
    and intrinsics from data table.
    Returns the constant table.
    """
    tsize = redshifts.shape[0]
    diag = numpy.square(_COEFF_SIGMAZ / redshifts)      # Peculiar redshift
    diag += numpy.square(redshifts * 0.055)             # Lensing
    for i in sm.range(tsize):                           # Intrinsic
        diag[i] += _SIGMA_COH_TABLE[setnum[i]] ** 2
    return diag


def loadcovbase(dirpath):
    """Read the "base" covariance matrix from the FITS files stored under
    the directory dirpath.
    """
    tmp = []
    for basename in os.listdir(dirpath):
        if basename.endswith(".fits"):
            fitsfile = fits.open(os.path.join(dirpath, basename))
            tmp.append(numpy.asarray(fitsfile[0].data))
            fitsfile.close()
    if not tmp:     # Data absent
        raise ValueError("cannot find sum of data covariance matrices")
    res = sum(tmp)
    return res


def loadsntable(path):
    """Read the supernova table from given data file located at path.
    Returns a numpy ndarray object, in which each row is a record
    and each column contain values for an attribute.  Only the useful ones
    are returned.
    """
    columns = (1,   # z, redshift
               4,   # m_b, peak B-magnitude
               6,   # X_1, shape parameter
               8,   # C, color parameter
               10,  # M_stellar, stellar mass of host galaxy
               17,  # dataset identifier
               )
    fulltable = numpy.loadtxt(path, usecols=columns)
    return fulltable


def _extend_order(order):
    """Extend an ordering of supernova objects to the corresponding ordering
    that indexes the (mag, stretch, color) triplet in the FITS matrices.
    """
    extended_order = []
    for i in order:
        t = i * 3
        extended_order.extend([t, t + 1, t + 2])
    return numpy.array(extended_order, dtype=int)


def _sort_table_and_basematrix(table, basematrix):
    """Sort the table by z and apply the order to basematrix.
    Returns sorted table and basematrix.
    """
    n = table.shape[0]
    order = numpy.argsort(table[:, 0], kind="mergesort")
    table_sorted = table[order]
    extended_order = _extend_order(order)
    basematrix_sorted = basematrix[numpy.ix_(extended_order, extended_order)]
    return table_sorted, basematrix_sorted


def _create_default_logbins():
    return numpy.asarray([numpy.linspace(-2, numpy.log10(1.3), 31)])


def _filter_items_by_bins(items, bins):
    """Filter the items by the given bins.

    Return value:
        A tuple of 2 parallel lists, where the 0th list is the indices of
        items, and the 1st is the corresponding binning result.  Those not
        belonging to any of the bins are excluded.
    """
    ids = []
    binresults = []
    for i, item in enumerate(items):
        binid = bins.searchenum(item)
        if binid[0] is not None:
            ids.append(i)
            binresults.append(binid)
    return (ids, binresults)


def _create_base_block(ids, basematrix):
    """Construct the "base block", a space-time trade-off for faster accessing
    the base matrix.

    Arguments:
        ids:  collection of data-point ids (integers).  Only the columns and
        rows corresponding to them are included in the output.

        basematrix:  the "base" matrix, as given by the binary data file
        themselves

    Return value:
        A numpy.ndarray object, B, of shape (3, 3, len(ids), len(ids)), where
        B[i][j] is the ij-th block matrix of covariances. i, j in {0, 1, 2},
        which corresponds to peak B-magnitude, stretch, and color parameters
        respectively.
    """
    mids = _extend_order(ids)
    section = basematrix[numpy.ix_(mids, mids)]
    bblock = numpy.empty((3, 3, len(ids), len(ids)))
    for i in sm.range(3):
        for j in sm.range(3):
            p = section[i::3, j::3]
            bblock[i][j] = p
    return numpy.ascontiguousarray(bblock)


class BinnedSN(object):
    """Binned supernova data.
    An instance of this class BinnedSN, after initialization, contains the
    following attributes that give access to the constants needed for
    computation:
        basematrix: Base covariance matrix for the observables,
                    corresponding to C_eta in B14.
        table: Array containing the data points.  Each row is a record
               for a point, and the columns are organized as follows:
                   0: z, redshift
                   1: m_b, peak B-magnitude
                   2: X_1, shape parameter
                   3: C, color parameter
                   4: M_stellar, stellar mass of host galaxy
                   5: dataset identifier
               Short-hand attributes are available for accessing copies
               of the table columns.  The columns has the following synonyms:
                   0: redshifts
                   1: peakbmags
                   2: shapes
                   3: colors
                   4: mstellar
                   5: dataset
        datadimension: Number of data points.
        bins: a binning.BinCollection instance corresponding to the binning
              scheme in log-redshift.
        binnings: List parallel to data columns in table, containing tuples of
                  the form (binid, (lower, upper)) for each corresponding data
                  point, where binid is the index of bin, and (lower, upper)
                  are the lower and upper ends of the bin.
        needhostcorr: Boolean array, precomputed for the predicate
                      M_stellar >= 10
        binidcontent: Mapping object that maps a bin id to a set
                      of data ids that are in the bin.
    """
    def __init__(self, basedirpath, tablepath, logbins=None, sort_by_z=False):
        """Initialize a BinnedSN instance by loading the data files
        and performing binning.

        basedirpath: Path to directory containing the FITS files
                     for various components in the base covariance
                     matrix.
        tablepath: Path to the text file containing the data table.
        logbins: Optional, ordered collection of ordered collection of numbers
                 that are the log-control points of the bins.  If not given,
                 the default binning method is used, which is the one adapted
                 to follow B14.
        sort: If True, the table and corresponding matrix will be sorted by
              redshift in ascending order (default: False)
        """
        if logbins is None:     # Use default binning scheme.
            self.bins = binning.BinCollection(_create_default_logbins())
        else:
            self.bins = binning.BinCollection(logbins)
        self._cpcache = list(self.bins.itercontrolpoints())
        basematrix = loadcovbase(basedirpath)
        # Symmetrize?
        basematrix = (basematrix + basematrix.T) / 2.0
        table = loadsntable(tablepath)
        if sort_by_z:
            table, basematrix = _sort_table_and_basematrix(table, basematrix)
        # Remove all data points that are out-of-bins.
        ids, binresults = _filter_items_by_bins(numpy.log10(table[:, 0]),
                                                self.bins)
        self.binnings = binresults
        self.table = numpy.ascontiguousarray(table[ids])
        # Corresponding ids in the matrix.
        self.datadimension = len(ids)
        if self.datadimension == 0:
            raise ValueError("No data selected.")
        self.bblock = _create_base_block(ids, basematrix)
        self.needhostcorr = self.table[:, 4] >= 10
        # Create "familiar" names, or views, for accessing table elements.
        self.redshifts = numpy.ascontiguousarray(self.table[:, 0])
        self.logredshifts = numpy.log10(self.redshifts)
        self.peakbmags = numpy.ascontiguousarray(self.table[:, 1])
        self.shapes = numpy.ascontiguousarray(self.table[:, 2])
        self.colors = numpy.ascontiguousarray(self.table[:, 3])
        self.mstellar = numpy.ascontiguousarray(self.table[:, 4])
        # This is not a view but a copy
        self.dataset = numpy.ascontiguousarray(self.table[:, 5], dtype=int)
        # Update the 00-block of bblock by the diagonal augment.
        # The meaning of this update is in Equation (13) of B14.
        self.bblock[0][0] += numpy.diag(_preeval_diag_in_cov(self.redshifts,
                                                             self.dataset))
        self.binidcontent = self._precompute_reverse_bin_lookup()
        return None

    def _precompute_reverse_bin_lookup(self):
        """Compute the reverse map of looking up the bin for given data id.

        Return value:
            A mapping of size nbins, where the key k maps to a collection of
            data ids that are in the kth bin.
        """
        rev = collections.defaultdict(set)
        for dataid in sm.range(self.datadimension):
            binid = self.binnings[dataid][0]
            rev[binid].add(dataid)
        return rev
