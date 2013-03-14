double *pyvector_to_Carraypxtrs(PyArrayObject *arrayin);
double **pymatrix_to_Carraypxtrs(PyArrayObject *arrayin);
double **ptrvector(long n);
void free_Carrayptrs(double **v);
static PyObject *build_Uab(PyObject *self, PyObject *args);
static PyObject *sum_dJ(PyObject *self, PyObject *args);