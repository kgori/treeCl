#!/usr/bin/env python

def extract_gamma_parameter(self):
        gamma_regex = re.compile(r'(?<=Gamma shape parameter: \t\t)[.\d+]+')
        try:
            gamma = float(gamma_regex.search(self.output).group())
        except AttributeError:
            print 'Couldn\'t extract parameters'
            return
        return gamma

def extract_GTR_parameters(self):
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
        Afreq = float(Afreq_regex.search(self.output).group())
        Cfreq = float(Cfreq_regex.search(self.output).group())
        Gfreq = float(Gfreq_regex.search(self.output).group())
        Tfreq = float(Tfreq_regex.search(self.output).group())
        AtoC = float(AtoC_regex.search(self.output).group())
        AtoG = float(AtoG_regex.search(self.output).group())
        AtoT = float(AtoT_regex.search(self.output).group())
        CtoG = float(CtoG_regex.search(self.output).group())
        CtoT = float(CtoT_regex.search(self.output).group())
        GtoT = float(GtoT_regex.search(self.output).group())
    except AttributeError:
        print 'Couldn\'t extract parameters'
        return

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