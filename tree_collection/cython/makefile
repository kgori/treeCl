all:	autowrap	install

autowrap:
	autowrap --out py_wrapper.pyx wrapper.pxd
	rm py_wrapper.cpp

build:
	python setup.py build_ext

inplace:
	python setup.py build_ext --inplace

install:
	python setup.py install

clean:
	rm -r build/ py_wrapper.* *.so