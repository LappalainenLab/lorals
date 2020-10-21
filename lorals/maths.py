#!/usr/bin/env python

import math
import time
import logging

from numbers import Real
from operator import itemgetter
from typing import Iterable, Iterator, Tuple, Union

def median(x: Iterable[Real]) -> Real:
    """Calculate the median
    >>> median([1, 2, 3])
    2
    >>> median([4, 3, 4, 1, 2])
    3
    """
    if isinstance(x, (int, float)):
        return x
    x_sorted: Tuple[Real, ...] = tuple(sorted(x))
    midpoint: int = int(len(x_sorted) / 2)
    if len(x_sorted) % 2:
        return x_sorted[midpoint]
    elif len(x_sorted) == 2:
        lower, upper = x_sorted # type: Real, Real
    else:
        lower, upper = x_sorted[midpoint:(midpoint + 2)] # type: Real, Real
    return upper - ((upper - lower) / 2)


def percentile(x: Iterable[Real], probs: Union[Real, Iterable[Real]]) -> Union[Real, Tuple[Real, ...]]:
    """Calculate a percentile
    >>> percentile(x = (15, 20, 35, 40, 50), probs = 40)
    29.0
    >>> percentile(x = (15, 20, 35, 40, 50), probs = (5, 40, 95))
    (16.0, 29.0, 48.0)
    """
    if isinstance(probs, (int, float)):
        probs: Tuple[Real] = (probs,)
    if not all(0 <= p <= 100 for p in probs):
        raise ValueError("All percentile values must be between 0 and 100, inclusive")
    xs: Tuple[Real, ...] = tuple(sorted(x))
    ranks: Tuple[Real, ...] = tuple(p / 100 * (len(xs) - 1) + 1 for p in probs)
    indexes: Tuple[int, ...] = tuple(int(math.floor(r)) for r in ranks)
    percs: Tuple[Real, ...] = tuple(xs[i - 1] + ((r % 1) * (xs[i] - xs[i - 1])) for r, i in zip(ranks, indexes))
    if len(percs) == 1:
        return percs[0]
    return percs


def pvalue_adjust(pvalues: Iterable[float]) -> Tuple[float, ...]:
    """Adjust p-values using the Benjamin-Hochberg method
    Basically, it's like R's p.adjust(method = 'fdr')
    >>> x = tuple(map(lambda x: round(x, 1), tuple(x * 0.1 for x in range(0, 10))))
    >>> pvalue_adjust(x)
    (0.0, 0.5, 0.6666666666666667, 0.75, 0.8, 0.8333333333333334, 0.8571428571428571, 0.875, 0.888888888888889, 0.9)
    """
    num_pvalues: int = len(pvalues)
    logging.info("Adjusting %s p-values", num_pvalues)
    adjust_start: float = time.time()
    order: Tuple[int, ...] = tuple(index for index, _ in sorted(enumerate(pvalues), key=itemgetter(1), reverse=True))
    ranked_order: itemgetter = itemgetter(*(index for index, _ in sorted(enumerate(order), key=itemgetter(1))))
    adj_value: Iterable[float] = (num_pvalues / index * pvalues[o] for index, o in zip(range(num_pvalues, 0, -1), order))
    adj_value: Tuple[float, ...] = tuple(min(value, 1.0) for value in adj_value)
    adj_value: Tuple[float, ...] = ranked_order(adj_value)
    logging.debug("Adjusting p-values took %s seconds", round(time.time() - adjust_start, 3))
    return adj_value
