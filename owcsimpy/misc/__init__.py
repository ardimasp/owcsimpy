from collections.abc import Iterable
import itertools

def flatten(arr):
    """ Flatten a messy iterable object into 
    an 1D iterable object.

    Examples
    --------
    >>> x = [[1,2],[1],[1,2,3,[4,5]]]
    >>> y = list(flatten(x))
    [1,2,1,1,2,3,4,5]

    
    """
    for i in arr:
        if isinstance(i,Iterable):
            yield from flatten(i)
        else:
            yield i