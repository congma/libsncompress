# vim: spell spelllang=en
"""
simplecache -- implementing cached methods using cachetools, with particular
interest in the methods that takes a numpy.ndarray argument.

This little module is built on top of cachetools
<http://pythonhosted.org/cachetools/>.  We provide an alternative to
cachetools.cachedmethod by the method decorator factory ``memoized()'', with
the particular target of a method that takes a numpy.ndarray argument.  It is
adapted from the implementation of cachetools.cachedmethod().  The original
cachetools package is developed by Thomas Kemmer (c) and is available as
MIT-Licensed free software, available from GitHub
<https://github.com/tkem/cachetools/> or PyPI.

Currently, it only supports methods with the definition signature of

    def method(self, arrayarg, *args, **kwargs):

where arrayarg is expected to be a numpy.ndarray instance, and it is this
argument that will be used to derive a key for cache access.  Class methods or
static methods are currently not supported.

For your class to benefit from this memoization decorator, it must first
inherit from our ArrayMethodCacheMixin class alongside with its other parent
classes.  After that, you can use the memoized() function to create decorators
that decorate your methods.  For example:
>>> import sys
>>> from six.moves import range
>>> import numpy
>>> class A(ArrayMethodCacheMixin, object):
...     v = 4.2
... 
...     @memoized()
...     def frob(self, array):
...         '''Docstring for method frob is preserved.'''
...         # Lengthy, expensive, and phony calculations...
...         tmp = array.copy()
...         p = numpy.outer(array, array)
...         for i in range(1000000):
...             tmp += numpy.dot(p, array)
...         return tmp + self.v
... 
...     @memoized(cachetype=cachetools.LFUCache, cachesize=2)
...     def spam(self, array, blah=5.0):
...         f = self.frob(array)
...         return numpy.dot(f, f + blah * self.v)
... 
...     @memoized(cachetype=None)
...     def eggs(self, array):
...         '''A special case for no-caching, for debug only.'''
...         return array * 2.0
... 
...     @memoized(cachetype=cachetools.LRUCache, cachesize=1,
...               getsizeof=sys.getsizeof)
...     def ham(self, array):
...         '''A very small cache.'''
...         print("ham: function body executed.")
...         return array

After that, you can interactively test the effect of memoization by
instantiating A:
>>> ta = A()
>>> x = numpy.array([1., 2., 3.])
>>> numpy.set_printoptions(precision=1)

The following call, when first invoked, will hang for a while:
>>> for v in ta.frob(x):
...     print(v)         # doctest: +NORMALIZE_WHITESPACE
14000005.2
28000006.2
42000007.2

But subsequent calls with the same argument will be very fast:
>>> for v in ta.frob(x):
...     print(v)         # doctest: +NORMALIZE_WHITESPACE
14000005.2
28000006.2
42000007.2

At this moment, if you wish, you can access the actual cache via the
_cachedict attribute, but manual handling of the cache is not recommended.
>>> ta._cachedict       # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
{'frob': LRUCache(..., maxsize=32, currsize=1)}

After the first call to ta.spam(), it will have its own cache, too:
>>> print("%.1f" % ta.spam(x))
2744002861600508.0
>>> items = list(ta._cachedict.items())
>>> items.sort()
>>> print(items)         # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
[('frob', LRUCache(..., maxsize=32, currsize=1)), ('spam', LFUCache(...))]

Docstring of the decorated method is preserved as-is:
>>> print(ta.frob.__doc__)
Docstring for method frob is preserved.
>>> print(ta.eggs.__doc__)
A special case for no-caching, for debug only.

If the cache is too small to store the array, it will be bypassed
transparently.
>>> _ = ta.ham(numpy.arange(4))
ham: function body executed.
>>> print(ta._cachedict["ham"])
LRUCache([], maxsize=1, currsize=0)
>>> _ = ta.ham(numpy.arange(4))
ham: function body executed.

Author: Cong Ma <cong.ma@obspm.fr>, (c) 2015.  See the file COPYING.
"""


import functools
import cachetools
import numpy


class NumKeyLite(object):
    """Read-Only ND-Array as Key (NumKeyLite).

    Should only be used as dict/set key and nothing else.
    """

    __slots__ = ("hashfcn", "__value", "__oldw", "__h", "__weakref__")

    def __init__(self, array, hashfcn=hash, *nparray_args, **nparray_kwargs):
        """NumKeyLite(array, *args, **kwargs) -- create key from
        numpy.ndarray instance other.

        The input array argument is converted to NumKeyLite instance.  To
        do this, numpy.asarray(array) is called first.  This may or may not
        copy the input argument.  See its docs for more.

        Extra positional and keyword arguments are passed to numpy.asarray().
        >>> numpy.set_printoptions(precision=1, floatmode="fixed")
        >>> a = [0, 1, 2]; ak = NumKeyLite(a, dtype=numpy.float64)
        >>> print(str(ak))  # doctest: +NORMALIZE_WHITESPACE
        [0.0 1.0 2.0]

        Side effect:  The underlying array is always set to read-only during
        the lifetime of the key instance.  When calling __init__(), if the
        "array" argument is not copied as a result of numpy.asarray(array)
        call, this array will be set to read-only by alias.  This is done to
        prevent data corruption.  User should fully understand the implication
        of using an ndarray as key, and should not attempt to manually set the
        array to writeable in this case.  Especially not during the creation of
        NumKeyLite objects!
        >>> from numpy.random import rand
        >>> a = rand(16); ak = NumKeyLite(a)
        >>> a[0] = 2.0
        Traceback (most recent call last):
            ...
        ValueError: assignment destination is read-only

        The writeable flag is restored after the key releases the reference:
        >>> del ak
        >>> a[0] = 2.0
        """
        self.hashfcn = hashfcn
        self.__value = numpy.asarray(array, *nparray_args, **nparray_kwargs)
        # As long as we're holding the reference to key, we make it read-only,
        # so that the key cannot be modified.
        self.__oldw = self.__value.flags.writeable
        self.__value.flags.writeable = False
        self.__h = hashfcn(self.__value.tobytes())
        return None

    def __hash__(self):
        return self.__h

    def __eq__(self, other):
        """Test whether the underlying data is bytewise-equal.
        >>> a = numpy.arange(16); b = a.copy()
        >>> ak = NumKeyLite(a); bk = NumKeyLite(b)
        >>> ak == bk
        True
        >>> ak != bk
        False
        >>> bk = NumKeyLite(b + 0.1)
        >>> ak == bk
        False
        >>> b = a.reshape(4, -1).copy(); bk = NumKeyLite(b)
        >>> ak == bk
        False
        >>> ak != bk
        True
        >>> ak == a  # doctest: +ELLIPSIS
        Traceback (most recent call last):
            ...
        TypeError: Not knowing how to compare.
        >>> ak = NumKeyLite(numpy.array(1.0))
        >>> bk = NumKeyLite(numpy.array(1.0))
        >>> ak == bk
        True
        """
        # NOTE: This is debatable.
        # NOTE: Hopefully, the hash function should behave so well that
        # collisions are rare, therefore making it rare that the expensive
        # equality comparison is called.
        try:
            # numpy.all() works on zero-dimensional arrays.
            return numpy.all(numpy.equal(self.__value, other.__value))
        except ValueError:
            return False
        except AttributeError:
            raise TypeError("Not knowing how to compare.")

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return str(self.__value)

    def __repr__(self):
        return "NumKeyLite(%r)" % self.__value

    def __del__(self):
        # When we (i.e. the NumKeyLite instance, not the value it holds) are
        # destroyed, restore old writeable flag to the referenced array.  This
        # may have nasty effects.
        self.__value.flags.writeable = self.__oldw


def memoized(cachetype=cachetools.LRUCache, cachesize=32,
             keyfcn=NumKeyLite, *cargs, **ckwargs):
    """Decorator-factory that returns a decorator suitable for methods of a
    class that inherits ArrayMethodCacheMixin.

    Optional arguments that fine-tunes the creation of method caches:
        cachetype: Type of cache -- I'd recommend the cache types from the
                   cachetools package.  It can actually be None, which means no
                   caching should be present.  This is probably only useful for
                   testing.
        cachesize: Size of cache.  The value of this argument will be the one
                   passed to the value of cachetype as the latter's first call
                   argument.  If the value of cachetype doesn't take a size
                   parameter this way, things may break.  This is the
                   convention followed by cachetools.  Except when cachetype is
                   None, then this argument is ignored.
        keyfcn: Key function to be applied to the array argument that is passed
                as the first mandatory argument of the method to be decorated.
        *cargs, **ckwargs: extra parameters to be passed to the value of
                           cachetype for instantiating the actual cache
                           instance.
    """
    if cachetype is None:
        # Do nothing, return identity decorator.
        return lambda x: x

    def decorator(arraymethod):
        """Method decorator to be returned by the enclosing factory function.
        """
        # Create the wrapper around the underlying method call.  This wrapper
        # does the memoization using the cache just
        # retrieve-if-absent-create'd.
        mname = arraymethod.__name__

        @functools.wraps(arraymethod)
        def wrapper(self, arrayarg, *args, **kwargs):
            """The actual wrapper that intercepts the method call arguments,
            performs the caching, and returns the result to caller.
            """
            # If no cache yet, create cache for this method.
            try:
                # _cachedict is keyed by the method names, rather than the
                # unwrapped method objects themselves, although the latter is
                # possible.  We choose the former because this helps debugging
                # better.  The original method, once wrapped, could be hard to
                # access by a Python name, although one can still enumerate
                # _cachedict's keys.  Our choice works, because
                # functools.wraps() ensures conservation of names.
                cache = self._cachedict[mname]
            except KeyError:
                cache = cachetype(cachesize, *cargs, **ckwargs)
                self._cachedict[mname] = cache
            argkey = keyfcn(arrayarg)
            try:
                return cache[argkey]
            except KeyError:
                # Cache miss, compute and store the return value.
                pass
            retval = arraymethod(self, arrayarg, *args, **kwargs)
            try:
                cache[argkey] = retval
            except ValueError:
                # Value probably too large, ignore and pass through.
                pass
            return retval
        return wrapper
    return decorator


class ArrayMethodCacheMixin(object):
    """Mix-in class that only creates an attribute to access the caches in the
    instance during initialization.

    Classes that inherits this mix-in can super()-delegate the __init__()
    method so that things happen automatically.
    """
    def __init__(self, *args, **kwargs):
        """Create the container for the method caches in this instance (self).
        """
        self._cachedict = {}
