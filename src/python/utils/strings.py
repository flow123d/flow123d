__author__ = 'Jan Hybs'



def join_iterable (iterable, prefix="", suffix="", separator=",", padding=None, extraSpace=2):
    result = ""
    result += prefix

    size = len (iterable)
    for i in range (size):
        if padding is None:
            result += str (iterable[i])
        else:
            result += str (iterable[i]).center(padding[i] + extraSpace)
        if i < size-1:
            result += (suffix+separator+prefix)
    result += suffix

    return result
