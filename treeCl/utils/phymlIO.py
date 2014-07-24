#!/usr/bin/env python
from __future__ import print_function
import re

"""
TO DO: Add PAML class to estimate these kinds of parameters better then
can dispense with this
"""


def extract_gamma_parameter(tree):
    gamma_regex = re.compile(r'(?<=Gamma shape parameter: \t\t)[.\d+]+')
    try:
        gamma = float(gamma_regex.search(tree.output).group())
    except AttributeError:
        print('Couldn\'t extract alpha parameter')
        return 1.0
    return gamma


def extract_GTR_parameters(tree):
    Afreq_regex = re.compile(r'(?<=f\(A\)= )[.\d+]+')
    Cfreq_regex = re.compile(r'(?<=f\(C\)= )[.\d+]+')
    Gfreq_regex = re.compile(r'(?<=f\(G\)= )[.\d+]+')
    Tfreq_regex = re.compile(r'(?<=f\(T\)= )[.\d+]+')
    AtoC_regex = re.compile(r'(?<=A <-> C    )[.\d+]+')
    AtoG_regex = re.compile(r'(?<=A <-> G    )[.\d+]+')
    AtoT_regex = re.compile(r'(?<=A <-> T    )[.\d+]+')
    CtoG_regex = re.compile(r'(?<=A <-> G    )[.\d+]+')
    CtoT_regex = re.compile(r'(?<=A <-> T    )[.\d+]+')
    GtoT_regex = re.compile(r'(?<=A <-> T    )[.\d+]+')

    try:
        Afreq = float(Afreq_regex.search(tree.output).group())
        Cfreq = float(Cfreq_regex.search(tree.output).group())
        Gfreq = float(Gfreq_regex.search(tree.output).group())
        Tfreq = float(Tfreq_regex.search(tree.output).group())
        AtoC = float(AtoC_regex.search(tree.output).group())
        AtoG = float(AtoG_regex.search(tree.output).group())
        AtoT = float(AtoT_regex.search(tree.output).group())
        CtoG = float(CtoG_regex.search(tree.output).group())
        CtoT = float(CtoT_regex.search(tree.output).group())
        GtoT = float(GtoT_regex.search(tree.output).group())
    except AttributeError:
        print('Couldn\'t extract GTR parameters')
        Afreq = Cfreq = Gfreq = Tfreq = 0.25
        AtoC = AtoG = AtoT = CtoG = CtoT = GtoT = 1.0

    d = dict(
        Afreq=Afreq,
        Cfreq=Cfreq,
        Gfreq=Gfreq,
        Tfreq=Tfreq,
        AtoC=AtoC,
        AtoG=AtoG,
        AtoT=AtoT,
        CtoG=CtoG,
        CtoT=CtoT,
        GtoT=GtoT,
    )

    return d
