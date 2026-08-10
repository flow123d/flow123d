#!/bin/python3
# author: Jan Hybs

import time


def ttl_cache(ttl=5.0):
    """
    Simple wrapper which caches last returned valued from the function for a given amount
    of time.
    Useful when dealing with expensive operations (such as system calls and more)

    Examples
    --------
    >>> import random
    ... random.seed(1234) # fix random generated numbers for clarification
    ... @ttl_cache(ttl=1.0)
    ... def foo():
    ...     return random.randint(1, 100)
    ...
    >>> foo() # function foo will be always called for the first time
    97
    >>> foo() # the following 1.0 second, foo won't be called
    97
    >>> time.sleep(1.0) # after that is will be called once again
    ... foo()
    ...
    45

    :param ttl: time to live (how long to treat the result as valid)
    :return:  wrapper
    """
    def ttl_decorator(func):
        func._last_called = 0
        func._result = None

        def func_wrapper(*args, **kwargs):
            diff = time.time() - func._last_called
            if diff > ttl:
                func._last_called = time.time()
                func._result = func(*args, **kwargs)
                return func._result
            else:
                return func._result
        return func_wrapper
    return ttl_decorator
