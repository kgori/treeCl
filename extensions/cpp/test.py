import nose
from nose import with_setup
import py_seq

def setup_func():
    global s
    s = py_seq.Seq('/Users/kgori/scratch/carnivores_16S.phy', 'phy', 
        'nt', True)

def teardown_func():
    pass

@with_setup(setup_func, teardown_func)
def test_load():
    assert isinstance(s, py_seq.Seq)

def test_copy():
    cop = s.__copy__()
    assert isinstance(cop, py_seq.Seq)

def test_bootstrap():
    boot = s.bootstrap_sample()
    assert isinstance(boot, py_seq.Seq)
    
def test_simulate():
    sim = s.simulate()
    assert isinstance(sim, py_seq.Seq)
    
