# vim: spell spelllang=en

import itertools
from operator import add
from bisect import bisect_left, bisect_right
import six.moves as sm


class BinCollection(object):
    """This is a collection of bins useful for binning the data into
    different windows.

    A collection is a list of bin-chains, in ascending order, in which each
    bin-chain does not intersect with any other.

    A bin-chain is an list of numbers (at least two), in ascending order, that
    defines a consecutive chain of bins.

    Should be slow but safe.

    The instance is created from a sequence of sequences of numbers.
    >>> a = [[1, 2, 3], [4, 5], [6, 7, 9, 12]]; c = BinCollection(a)

    After initialization, information about bin numbers are saved.
    >>> c.nchains
    3
    >>> c.nbins
    6

    To search for the bin to which a number belongs, use searchenum(),
    which finds the leftmost bin to which the number resides in.
    >>> c.searchenum(2.2)
    (1, (2, 3))
    >>> c.searchenum(4.2)
    (2, (4, 5))
    >>> c.searchenum(9.00001)
    (5, (9, 12))
    >>> c.searchenum(1)
    (0, (1, 2))
    >>> c.searchenum(4)
    (2, (4, 5))
    >>> c.searchenum(5)
    (2, (4, 5))
    >>> c.searchenum(6.0)
    (3, (6, 7))
    >>> c.searchenum(7)
    (3, (6, 7))
    >>> c.searchenum(0.99999) == (None, (None, None))
    True
    >>> c.searchenum(5.5) == (None, (None, None))
    True
    >>> c.searchenum(12.1) == (None, (None, None))
    True

    To split a bin index into the chain and bin-in-chain parts, use
    binaddress()
    >>> c.binaddress(0)
    (0, 0)
    >>> c.binaddress(2)
    (1, 0)
    >>> c.binaddress(4)
    (2, 1)
    >>> c.binaddress(5)
    (2, 2)
    >>> c.binaddress(2.5)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: Invalid index value: ...

    To get a bin by its index, use getbin()
    >>> c.getbin(0)
    (1, 2)
    >>> c.getbin(5)
    (9, 12)
    >>> c.getbin(-5)
    (2, 3)
    >>> c.getbin(100)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    IndexError: Index out of range.
    >>> c.getbin(-7)  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    IndexError: Index out of range.

    To iterate over control points, use itercontrolpoints()
    >>> [i for i in c.itercontrolpoints()]
    [1, 2, 3, 4, 5, 6, 7, 9, 12]
    >>> [i for i in c.itercontrolpoints_classify()]  # doctest: +ELLIPSIS
    [(1, None, 0), (2, 0, 1), (3, 1, None), (4, None, 2), (5, 2, None)...]

    >>> c.getbinendids(1)
    (1, 2)
    >>> c.getbinendids(4)
    (6, 7)
    """
    def __init__(self, otherlist):
        """Initialize using the list otherlist.

        Validation is performed to ensure non-overlapping bin-chains.
        Each item in otherlist is itself a list that defines a bin-chain.
        If this item is not in ascending order, it is sorted.

        ValueError is raised if any item in otherlist cannot be used as
        a valid bin-chain (e.g. insufficient number of items, zero-length
        interval, etc.)
        >>> c = BinCollection([[1, 2], [4]])  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        ValueError: Chain too small: ...
        >>> c = BinCollection([[1, 2, 2, 3], [4]])  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        ValueError: Duplicate item in chain: ...

        ValueError is also raised if the collection cannot be constructed
        due to overlapping bin-chains.
        >>> c = BinCollection([[1, 2], [4, 5], [4.5, 6]])  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        ValueError: Overlapping chains in input: ...
        """
        vchains = [validatechain(chain) for chain in otherlist]
        # 1st pass: sort in ascending order
        vchains.sort(key=lambda l: l[0])
        # 2nd pass: checking for overlapping.  We even reject the chains
        # that could have been merged.  Not worth the trouble now ;)
        # Essentially we're asking for strict non-overlapping.
        for left, right in pairwise(vchains):
            if areoverlapping(left, right):
                raise ValueError("Overlapping chains in input: %s" % vchains)
        self.store = tuple(vchains)
        self.nchains = len(vchains)
        self.chainlens = [len(chain) - 1 for chain in vchains]
        self.nbins = sum(self.chainlens)
        self.ncontrolpoints = self.nbins + self.nchains
        b_cache = list(self.iterbins())
        self._lefts_cache = [x[0] for x in b_cache]
        self._rights_cache = [x[1] for x in b_cache]
        return None

    def itercontrolpoints(self):
        """Returns iterator that iterates through all control points."""
        return itertools.chain.from_iterable(self.store)

    def itercontrolpoints_classify(self):
        """Returns iterator that yields the control points along with
        their status of being inner nodes or chain-end nodes.
        Each element yielded is a three-tuple (point, leftbin, rightbin).
        point is the value of the control point.
        leftbin and rightbin are the bin id of the bin to the left/right
        of this control point.  They can be None, indicating that the point
        is a left/right chain-end point.
        """
        for chainid, chain in enumerate(self.store):
            for pointid, point in enumerate(chain):
                right = sum(self.chainlens[:chainid]) + pointid
                left = right - 1
                if pointid == 0:
                    left = None
                if pointid == self.chainlens[chainid]:
                    right = None
                yield (point, left, right)

    def iterbins(self):
        """Iterator that iterates through all bins (given by 2-tuples).
        """
        for chain in self.store:
            for pair in pairwise(chain):
                yield pair

    def searchenum(self, item):
        """Search for the location of item in the bins.  If found, return
        (binindex, (bin_lower, bin_upper));  if not, return
        (None, (None, None)).
        """
        # Cf: Python Library Reference, "bisect", "Searching Sorted Lists"
        i = bisect_left(self._rights_cache, item)
        if i != self.nbins and item >= self._lefts_cache[i]:
            return (i, (self._lefts_cache[i], self._rights_cache[i]))
        return (None, (None, None))

    def binaddress(self, index):
        """Get the "address" for the given bin index.
        The address is a tuple (chainid, relbinid),
        where chainid is the index of the chain, and relbinid
        is the relative id of this bin in the chain.
        """
        try:
            cl_cache = self._cl_cache
        except AttributeError:
            cl_cache = list(scanl1(add, self.chainlens))
            self._cl_cache = cl_cache
        index = self._normalize_bin_index(index)
        i = bisect_right(cl_cache, index)
        relid = self.chainlens[i] - (cl_cache[i] - index)
        return i, relid

    def getbin(self, index):
        """Get the nth bin given index."""
        chainid, relbinid = self.binaddress(index)
        chain = self.store[chainid]
        return (chain[relbinid], chain[relbinid + 1])

    def getbinendids(self, index):
        """Get the indices of the bin's end-points, viewed as member
        of a flattened list of control points, given bin index.
        """
        chainid, relbinid = self.binaddress(index)
        # Can be proved using induction
        lower = sum(self.chainlens[0:chainid]) + chainid + relbinid
        return (lower, lower + 1)

    def _normalize_bin_index(self, index):
        """Return normalized bin index.
        """
        if index != int(index):
            raise ValueError("Invalid index value: %r" % index)
        if index >= self.nbins or index <= -1 - self.nbins:
            raise IndexError("Index out of range.")
        if index < 0:
            index += self.nbins
        return index


def validatechain(chaincandidate):
    """Validate a single chain, make it sorted.  Raise ValueError
    for an invalid chain.
    """
    sortedchain = sorted(chaincandidate)
    lchain = len(sortedchain)
    # Check for minimal chain length
    if lchain < 2:
        raise ValueError("Chain too small: %s" % sortedchain)
    # Check for duplication within chain
    for i in sm.range(lchain - 1):
        if sortedchain[i] == sortedchain[i + 1]:
            raise ValueError("Duplicate item in chain: %s" %
                             sortedchain)
    return sortedchain


def pairwise(iterable):
    """Taken straight out of Python Library Reference.
    s -> (s0, s1), (s1, s2), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return sm.zip(a, b)


def scanl1(func, seq):
    """Shamelessly copied from Haskell.
    >>> a = sm.range(5)
    >>> list(scanl1(add, a))
    [0, 1, 3, 6, 10]
    """
    it = iter(seq)
    res = next(it)
    yield res
    for thing in it:
        res = func(res, thing)
        yield res


def areoverlapping(left, right):
    """Test if the chains represented by left and right are overlapping."""
    return right[0] <= left[-1]
