#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include "evrot_extensions.h"

// Utilities
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin) {
    int n;

    n=arrayin->dimensions[0];
    return (double *) arrayin->data;
}

double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin) { 
   double **c, *a;
   int i,n,m;

   n=arrayin->dimensions[0];
   m=arrayin->dimensions[1];
   c=ptrvector(n);
   a=(double *) arrayin->data; /* pointer to arrayin data as double */
   for ( i=0; i<n; i++) {
      c[i]=a+i*m; }
   return c;
}

double **ptrvector(long n) { 
   double **v;
   v=(double **)malloc((size_t) (n*sizeof(double)));
   if (!v)   {
      printf("In **ptrvector. Allocation of memory for double array failed.");
      exit(0); }
   return v;
}

void free_Carrayptrs(double **v) { 
   free((char*) v);
}

// Exported methods
static PyObject *sum_dJ(PyObject *self, PyObject *args) {
    // arguments
    PyArrayObject *A_x_Y, *Y_sq;                 // 2D arrays
    PyArrayObject *mv_sq, *mv_cb, *max_A_values; // 1D arrays
    int dim, ndata;

    // locals
    int i, j;                                   // array indexing
    double **p_A_x_Y, **p_Y_sq;                 // pointers to 2D arrays
    double *p_mv_sq, *p_mv_cb, *p_max_A_values; // pointers to 1D arrays
    double tmp1, tmp2, dJ;

    // parse input
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!ii",
        &PyArray_Type, &A_x_Y,
        &PyArray_Type, &Y_sq,
        &PyArray_Type, &mv_sq,
        &PyArray_Type, &mv_cb,
        &PyArray_Type, &max_A_values,
        &dim, &ndata)) return NULL; 

    // type checking
    if (A_x_Y->nd != 2 || A_x_Y->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "A_x_Y must be a two-dimensional float array");
        return NULL;
    }
    if (Y_sq->nd != 2 || Y_sq->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "Y_sq must be a two-dimensional float array");
        return NULL;
    }
    if (mv_cb->nd != 1 || mv_cb->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "mv_cb must be a one-dimensional float array");
        return NULL;
    }
    if (mv_sq->nd != 1 || mv_sq->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "mv_sq must be a one-dimensional float array");
        return NULL;
    }
    if (max_A_values->nd != 1 || max_A_values->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "max_A_values must be a one-dimensional float array");
        return NULL;
    }

    // set pointers
    p_A_x_Y = pymatrix_to_Carrayptrs(A_x_Y); // remember to free this memory!
    p_Y_sq = pymatrix_to_Carrayptrs(Y_sq);   // ditto
    p_mv_sq = pyvector_to_Carrayptrs(mv_sq);
    p_mv_cb = pyvector_to_Carrayptrs(mv_cb);
    p_max_A_values = pyvector_to_Carrayptrs(max_A_values);

    dJ = 0;
    for (j = 0; j < dim; j++) {
        for (i = 0; i < ndata; i++) {
            tmp1 = (p_A_x_Y[i][j]) / (p_mv_sq[i]);
            tmp2 = (p_max_A_values[i]) * (p_Y_sq[i][j]) / (p_mv_cb[i]);
            dJ += tmp1 - tmp2;
        }
    }
    dJ = 2*dJ/ndata/dim;

    // free allocated memory from 2D array pointers
    free_Carrayptrs(p_A_x_Y);
    free_Carrayptrs(p_Y_sq);
    
    return Py_BuildValue("d",dJ);

}

static PyObject *build_Uab(PyObject *self, PyObject *args) {
    
    // arguments
    PyArrayObject *theta, *ik, *jk;
    int a, b, dim;

    // locals
    PyArrayObject *Uab;                     // Output 2D array
    npy_intp i, j, k, m1, m2, n;            // np array indexing
    npy_long *p_m1, *p_m2;                  // pointers to integer values
    double *p_theta, **p_Uab; // pointers to float values
    double tmp, cos_theta, sin_theta;       // float data
    // int dims[2];
    npy_intp dims[2];

    // parse input
    if (!PyArg_ParseTuple(args, "O!iiO!O!i",
        &PyArray_Type, &theta, 
        &a, &b, 
        &PyArray_Type, &ik,
        &PyArray_Type, &jk,
        &dim)) return NULL;

    // type check on input
    if (theta->nd != 1 || theta->descr->type_num != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "Theta must be a one-dimensional float array");
        return NULL;
    }
    if (ik->nd != 1 || ik->descr->type_num != NPY_LONG) {
        PyErr_SetString(PyExc_ValueError,
        "ik must be a one-dimensional integer array");
        return NULL;
    }
    if (jk->nd != 1 || jk->descr->type_num != NPY_LONG) {
        PyErr_SetString(PyExc_ValueError,
        "jk must be a one-dimensional integer array");
        return NULL;
    }

    dims[0] = dims[1] = dim;
    // Uab = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE); // 2D Array
    Uab = (PyArrayObject *) PyArray_SimpleNew(2,dims,NPY_DOUBLE); // 2D Array

    // pointers
    p_theta = pyvector_to_Carrayptrs(theta);
    p_Uab = pymatrix_to_Carrayptrs(Uab); // remember to free this memory!
    
    // make Identity matrix of size dim
    for (i=0; i < dim; i++){
        // set diagonal to 1
        p_Uab[i][i] = 1.0;
        // set off-diagonal to 0
        for (j=i+1; j < dim; j++ ){
            p_Uab[i][j] = 0.0;
            p_Uab[j][i] = 0.0;
        }
    }

    if (b < a) {

        free_Carrayptrs(p_Uab); // free 2D array pointer
        return PyArray_Return(Uab);
    }

    for (k=a; k<=b; k++){
        // values that change with k
        p_m1 = (npy_long *)(ik->data + k*ik->strides[0]); // pointer to ik[k]
        p_m2 = (npy_long *)(jk->data + k*jk->strides[0]); // pointer to jk[k]
        m1 = *p_m1; // the value at ik[k]
        m2 = *p_m2; // the value at jk[k]
        cos_theta = cos(p_theta[k]);
        sin_theta = sin(p_theta[k]);
        for (n=0; n<dim; n++){
            // values that change with n
            tmp = p_Uab[n][m1] * cos_theta - p_Uab[n][m2] * sin_theta;
            p_Uab[n][m2] = p_Uab[n][m1] * sin_theta + p_Uab[n][m2] * cos_theta;
            p_Uab[n][m1] = tmp;
        }
    }
    free_Carrayptrs(p_Uab); // free 2D array pointer
    return PyArray_Return(Uab); 
}

static PyMethodDef evrot_extensions_methods[] = {
    {"build_Uab", build_Uab, METH_VARARGS},
    {"sum_dJ", sum_dJ, METH_VARARGS}
};

void initevrot_extensions() {
    (void) Py_InitModule("evrot_extensions", evrot_extensions_methods);
    import_array();
}
