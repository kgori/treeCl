from __future__ import print_function

# standard library
from copy import copy
from subprocess import Popen, PIPE

try:
    from subprocess import DEVNULL
except ImportError:
    import os

    DEVNULL = open(os.devnull, 'wb')

# third party
import numpy as np

# treeCl
from external import ExternalSoftware


def guess_seqtype(rec, threshold=0.8):
    """ Looks at proportion of As, Cs, Gs and Ts across all
    sequences in the record, and guesses the data to be dna
    if the proportion is > threshold, else it guesses protein """

    newrec = copy(rec)
    newrec.change_case('upper')
    all_seqs = ''.join(newrec.sequences)
    nuc_freqs = sum(all_seqs.count(x) for x in ['A', 'C', 'G', 'T'])
    nuc_freqs /= float(len(all_seqs))
    return 'dna' if nuc_freqs > threshold else 'protein'


class DV(ExternalSoftware):
    default_binary = 'darwin'

    def __init__(self, record, verbosity=0):
        super(DV, self).__init__(tmpdir='')
        self.record = record
        if not self.record.datatype:
            self.record.datatype = guess_seqtype(self.record)
        self.fasta = "[ ['{}'], ['{}'] ]".format('\',\''.join(record.sequences),
                                                 '\',\''.join(record.headers))
        self.cmd = '''# Darwin Script to make TreeCollection work
# Need to create an MAlignment object from sequence alignment file
# Then call RobustEstimateDistVarM on the MAlignment object,
# RobustEstimateDistVarM calls EstimatePamNoGap, which in turn calls RemoveGaps
# JUST RETURNS DISTANCE-VARIANCE MATRIX

# TO BRING IN COMMAND LINE ARGS TO DARWIN USE THIS CALL:
# SHELL>$ echo "fil := ReadFastaWithNames('<infile.fas>'); seqtype := 'AA'/'DNA'; fpath := '<output_filepath>'; ReadProgram('<path_to_this_file>/TC_wrapper.drw'); " | darwin


CreateDNAMatrices := proc( cnts:matrix(numeric) ;
    'FD0' = ((FD0=-28.045838538909323021):numeric),
    'FD1' = ((FD1=10.347906835124965854):numeric),
    'ID0' = ((ID0=-.43546673485399587638):numeric) )
global DMS;

#                      [.00050332799790683216930     -.00035779222938728085945]
#  CovarianceFD0FD1 := [                                                      ]
#                      [-.00035779222938728085945    .00025554907817687981940 ]
#                                                           -7
#                    VarianceID0 := .25143247825684950187 10
#
#                   StdDeviationID0 := .00015856622536241742230

if nargs = 0 then
     CreateDayMatrices(
    [[25212428,0,0,0,1072392,0,0,3666466,0,0,0,0,0,0,0,0,977696.5,0,0,0],
     [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [1072392,0,0,0,20273981,0,0,909400.5,0,0,0,0,0,0,0,0,3672543,0,0,0],
     [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
     [3666466,0,0,0,909400.5,0,0,21654018,0,0,0,0,0,0,0,0,1062248,0,0,0],
     [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
     [977696.5,0,0,0,3672543,0,0,1062248,0,0,0,0,0,0,0,0,27868216,0,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]
 )
else CreateDayMatrices(cnts) fi;

if FD0 >= 0 or FD0+FD1 >= 0 then
    error('FD0 and FD1  will cause deletion costs to be positive') fi;
DMS2 := [];
p := 1/10000;
while p < 5 do DMS2 := append(DMS2,DayMatrix(p));  p := p*1.02 od:
for p from ceil(10*DMS2[-1,PamNumber]) to 1500 do
    DMS2 := append(DMS2,DayMatrix(p/10))
od:
for i to length(DMS) while DMS[i,PamNumber] <= 150 do od;
DMS2 := [op(DMS2),op(DMS[i..-1])];

for z in DMS2 do
    z[FixedDel] := min(-1,FD0 + FD1*log10(z[PamNumber]));
    z[IncDel] := ID0;
od:

DMS := DMS2;
NULL;
end:

RemoveGaps := proc(t:string, u:string)
    if length(t) <> length(u) then error('inconsistent lengths'); fi;
    tnew := [];
    unew := [];
    for i to length(t) do
        if t[i] <> '_' and u[i] <> '_' then
           tnew := append(tnew,t[i]);
           unew := append(unew,u[i]);
        fi;
    od;
    return(string(op(tnew)),string(op(unew)));
end:

EstimatePamNoGap := proc(s1,s2,myDMS)
    if s1 = '' or s2 = '' then
        return([0,100,0]);
    fi;
    tmp := [RemoveGaps(s1,s2)];
    if tmp[1] = '' or tmp[2] = '' then
        return([0,100,0]);
    fi;
    res := traperror(EstimatePam(op(tmp),myDMS));
    if res = lasterror then
        return([0,100,0]);
    else
        return(res);
    fi;
end:

RobustEstimateDistVarM :=  proc(dmsa:MAlignment, seqtype:{'AA','DNA'})
    global DMS;

    if not assigned(DMS) then
        if seqtype = 'AA' then
            CreateDayMatrices();
        else
            CreateDNAMatrices();
            # CreateDayMatrices(HumanMtDNA);
        fi;
    else
        if seqtype = 'AA' and DMS[1,Sim,2,2] = -50 then
            CreateDayMatrices();
        elif seqtype = 'DNA' and DMS[1,Sim,2,2] <> -50 then
            CreateDNAMatrices();
            # CreateDayMatrices(HumanMtDNA);
        fi;
    fi;
    msa := dmsa['AlignedSeqs'];
    ls := length(msa);
    for k to ls do
        msa[k] := ReplaceString('-','_',msa[k]);
        if seqtype = 'DNA' then
            msa[k] := ReplaceString('N','X',msa[k]);
        fi;
    od;

    Dmsa := CreateArray(1..ls, 1..ls, 0);
    Vmsa := CreateArray(1..ls, 1..ls, 0);
    for k to ls do for l from k+1 to ls do
        tmp := EstimatePamNoGap(msa[k], msa[l], DMS);
        Dmsa[k,l] := Dmsa[l,k] := tmp[2] / 100.0; # PAM -> subs / site
        Vmsa[k,l] := Vmsa[l,k] := tmp[3] / 10000.0; # PAM -> subs / site
    od od;

    return([Dmsa,Vmsa]);
end:

Tri := proc(mtx1:array, mtx2:array)
    assert(length(mtx1)=length(mtx2));
    n := length(mtx1):
    #ret = CreateArray( 1..n*(n-1)*0.5 ):
    ret := [];
    for i to n do
        row := [];
        for j from 1 to i do
            row := append(row,mtx2[i][j]);
        od;
        for k from i+1 to n do
            row := append(row,mtx1[i][k]);
        od;
        ret := append(ret,row);
    od;
    return( ret );
end:

TrimLabels := proc(labels:array)
    assert(length(labels) > 0);
    len := length(labels);
    ret := [];
    for i to len do
        newlab := [];
        for j to length(labels[i]) do
            if labels[i][j] <> '/' then
                newlab := append(newlab, labels[i][j]);
            else
                break;
            fi
        od;
        newlab := string(op(newlab));
        ret := append(ret, newlab);
    od;
    return( ret );
end:

####################################################################################################

fil := %s;
#print(fil);
alignedSeqs := CreateArray(1..length(fil[1])):
labs := TrimLabels(fil[2]):
ntaxa := length(labs):
#print(labs);
unalignedSeqs := CreateArray(1..length(fil[1])):
for i to length(fil[1]) do
    alignedSeqs[i] := uppercase(ReplaceString('-','_',fil[1][i])):
    unalignedSeqs[i] := uppercase(ReplaceString('-','',fil[1][i])):
od:
MSA := MAlignment(unalignedSeqs,alignedSeqs,labs):
dvm := RobustEstimateDistVarM(MSA,%s):
tri := Tri(dvm[1],dvm[2]):
# PRINT DistVar.txt ###########################################################
#print ('DistVar.txt');
dv := ConcatStrings( tri[1], ' ');
for i from 2 to length(tri) do
    line := ConcatStrings( tri[i], ' ');
    dv := ConcatStrings( [ dv, line ], '\n');
od;
WriteData(dv, terminal);
quit;
''' % (self.fasta, ('AA' if self.record.datatype == 'protein' else 'DNA'))

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, record):
        newrec = copy(record)
        newrec.change_case('upper')
        self._record = newrec

    def execute(self):
        p = Popen(self.binary, stdout=PIPE, stdin=PIPE, stderr=DEVNULL)
        return p.communicate(input=self.cmd)[0]

    def parse(self, output):
        lines = output.rstrip().split('\n')
        end = lines.index('> WriteData(dv, terminal);')
        lines = lines[end + 1:]
        return np.array([line.split() for line in lines], dtype=float)

    def run(self):
        s = self.execute()
        return self.parse(s)


def runDV(record, verbosity=0):
    dw = DV(record, verbosity)
    dv_matrix = dw.run()
    labels = ' '.join(record.headers)
    record.dv.append((dv_matrix, labels))
    return (dv_matrix, labels)
