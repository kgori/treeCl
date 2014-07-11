all: autowrap setup test

autowrap: seq.pxd
	autowrap --out py_seq.pyx seq.pxd

setup:
	CC=gcc CXX=g++ python setup.py build_ext --inplace
	# /usr/bin/python setup.py build_ext --inplace

test:
	nosetests test.py

clean:
	rm -r build/ py_seq.*
