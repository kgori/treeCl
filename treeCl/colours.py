from __future__ import division
from builtins import range
import math
import collections

# http://stackoverflow.com/a/20715034

def hcl_to_rgb(hue=0, chroma=0, luma=0) :
    # Notes:
    #   coded from http://en.wikipedia.org/wiki/HSL_and_HSV#From_luma.2Fchroma.2Fhue
    #   with insights from gem.c in MagickCore 6.7.8
    #   http://www.imagemagick.org/api/MagickCore/gem_8c_source.html
    # Assume:
    #   h, c, l all in range 0 .. 1 (cylindrical coordinates)
    # Returns a tuple:
    #   r, g, b all in the range 0 .. 1 (cubic cartesian coordinates)

    # sanity checks
    hue = math.modf(float(hue))[0]
    if hue < 0 or hue >= 1 :
        raise ValueError('hue is a value greater than or equal to 0 and less than 1')
    chroma = float(chroma)
    if chroma < 0 or chroma > 1 :
        raise ValueError('chroma is a value between 0 and 1')
    luma = float(luma)
    if luma < 0 or luma > 1 :
        raise ValueError('luma is a value between 0 and 1')

    # do the conversion
    _h = hue * 6.0
    x = chroma * ( 1 - abs((_h % 2) - 1) )

    c = chroma
    if   0 <= _h and _h < 1 :
        r, g, b = (c, x, 0.0)
    elif 1 <= _h and _h < 2 :
        r, g, b = (x, c, 0.0)
    elif 2 <= _h and _h < 3 :
        r, g, b = (0.0, c, x)
    elif 3 <= _h and _h < 4 :
        r, g, b = (0.0, x, c)
    elif 4 <= _h and _h < 5 :
        r, g, b = (x, 0.0, c)
    elif 5 <= _h and _h <= 6 :
        r, g, b = (c, 0.0, x)
    else :
        r, g, b = (0.0, 0.0, 0.0)

    m = luma - (0.298839*r + 0.586811*g + 0.114350*b)
    z = 1.0
    if m < 0.0 :
        z = luma/(luma-m)
        m = 0.0
    elif m + c > 1.0 :
        z = (1.0-luma)/(m+c-luma)
        m = 1.0 - z * c
    (r, g, b) = (z*r+m, z*g+m, z*b+m)

    # clipping ...
    (r, g, b) = (min(r, 1.0), min(g, 1.0), min(b, 1.0))
    (r, g, b) = (max(r, 0.0), max(g, 0.0), max(b, 0.0))
    return (r, g, b)

def ggColorSlice(n=12, hue=(0.004,1.00399), chroma=0.8, luma=0.6, skipHue=True) :
    # Assume:
    #   n: integer >= 1
    #   hue[from, to]: all floats - red = 0; green = 0.33333 (or -0.66667) ; blue = 0.66667 (or -0.33333)
    #   chroma[from, to]: floats all in range 0 .. 1
    #   luma[from, to]: floats all in range 0 .. 1
    # Returns a list of #rgb colour strings:

    # convert stand alone values to ranges
    if not isinstance(hue, collections.Iterable):
        hue = (hue, hue)
    if not isinstance(chroma, collections.Iterable):
        chroma = (chroma, chroma)
    if not isinstance(luma, collections.Iterable):
        luma = (luma, luma)

    # convert ints to floats
    hue = [float(hue[y]) for y in (0, 1)]
    chroma = [float(chroma[y]) for y in (0, 1)]
    luma = [float(luma[y]) for y in (0, 1)]

    # some sanity checks
    n = int(n)
    if n < 1 or n > 360 :
        raise ValueError('n is a value between 1 and 360')
    if any([chroma[y] < 0.0 or chroma[y] > 1.0 for y in (0, 1)]) :
        raise ValueError('chroma is a value between 0 and 1')
    if any([luma[y] < 0.0 or luma[y] > 1.0 for y in (0, 1)]) :
        raise ValueError('luma is a value between 0 and 1')

    # generate a list of hex colour strings
    x = n + 1 if n % 2 else n
    if n > 1 :
        lDiff = (luma[1] - luma[0]) / float(n - 1.0)
        cDiff = (chroma[1] - chroma[0]) / float(n - 1.0)
        if skipHue :
            hDiff = (hue[1] - hue[0]) / float(x)
        else :
            hDiff = (hue[1] - hue[0]) / float(x - 1.0)
    else:
        hDiff = 0.0
        lDiff = 0.0
        cDiff = 0.0

    listOfColours = []
    for i in range(n) :
        c = chroma[0] + i * cDiff
        l = luma[0] + i * lDiff
        h = math.modf(hue[0] + i * hDiff)[0]
        h = h + 1 if h < 0.0 else h
        (h, c, l) = (min(h, 0.99999999999), min(c, 1.0), min(l, 1.0))
        (h, c, l) = (max(h, 0.0), max(c, 0.0), max(l, 0.0))
        (r, g, b) = hcl_to_rgb(h, c, l)
        listOfColours.append( '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) )
    return listOfColours
