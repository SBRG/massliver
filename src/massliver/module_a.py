def a():
    """An example function 'a'

    Returns
    -------
    str
        The letter 'a'.
    """
    return "a"


def b():
    """An example function 'b'

    Returns
    -------
    bool
    """
    if a() == "a":
        return True
    else:
        print("Not true")
    return False


__all__ = [
    "a",
    "b",
]
