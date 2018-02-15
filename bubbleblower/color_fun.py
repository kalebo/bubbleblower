PHI_CONJ = 0.618033988749895

# Derived from Martin Ankerl's blog about random programmatic colors.
def int_to_hsv(i: int, to_string=False):
    """Maps an integer to a color in hvs. Colors are confirmed to be unique in the
    range [1, 100] however the optimal visually distinguishable colors lay in
    periods of about 16
    """

    h = ((i * .944) + PHI_CONJ) % 1
    s = ((i * .7) + PHI_CONJ) % 1
    v = ((i * .1) + PHI_CONJ) % 1
    sa = s if (s > .5) else s + .4
    va = v if (v > .5) else v + .4

    hsv = (h, sa, va)
    if to_string:
        return " ".join(map(str, hsv))
    return hsv
