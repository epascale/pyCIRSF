#include <stdio.h>
#include <math.h>
#include "Python.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "numpy/arrayobject.h"

static PyObject *medfilt2d(PyObject *self, PyObject *args);

static PyMethodDef _medfilt[] = {
	{"medfilt2d", medfilt2d, METH_VARARGS, "compute median filter"},
        {NULL, NULL}
};




int  not_doublematrix(PyArrayObject *mat)  {
	if (PyArray_TYPE(mat) != NPY_DOUBLE || PyArray_NDIM(mat) != 2)  {
		PyErr_SetString(PyExc_ValueError,
			"In not_doublematrix: array must be of type Float and 2 dimensional (n x m).");
		return 1;  }
	return 0;
}

int  not_boolmatrix(PyArrayObject *mat)  {
	if (PyArray_TYPE(mat) != NPY_BOOL || PyArray_NDIM(mat) != 2)  {
		PyErr_SetString(PyExc_ValueError,
			"In not_boolmatrix: array must be of type Bool and 2 dimensional (n x m).");
		return 1;  }
	return 0;
}

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double nr_select(int k, int n, double arr[])
{
    int i,ir,j,l,mid;
    double a,temp;

    l=1;
    ir=n;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[ir] < arr[l]) {
                    SWAP(arr[l],arr[ir])
            }
            return arr[k];
        } else {
            mid=(l+ir) >> 1;
            SWAP(arr[mid],arr[l+1])
            if (arr[l] > arr[ir]) {
                    SWAP(arr[l],arr[ir])
            }
            if (arr[l+1] > arr[ir]) {
                    SWAP(arr[l+1],arr[ir])
            }
            if (arr[l] > arr[l+1]) {
                    SWAP(arr[l],arr[l+1])
            }
            i=l+1;
            j=ir;
            a=arr[l+1];
            for (;;) {
                    do i++; while (arr[i] < a);
                    do j--; while (arr[j] > a);
                    if (j < i) break;
                    SWAP(arr[i],arr[j])
            }
            arr[l+1]=arr[j];
            arr[j]=a;
            if (j >= k) ir=j-1;
            if (j <= k) l=i;
        }
    }
}
#undef SWAP

static PyObject *medfilt2d(PyObject *self, PyObject *args)
{
    PyArrayObject *matin, *maskin, *matout;
    double *cin, *cout, *buf;
    unsigned char *mask;
    int wy, wx, offs_x, offs_y, step;
    int ny, nx, dims[2];
    int i, j, x, y, n, idx;
    
    /* Parse tuples separately since args will differ between C fcns */
    
    if (!PyArg_ParseTuple(args, "O!O!iii", 
            &PyArray_Type, &matin, &PyArray_Type, &maskin, &wy, &wx, &step))  
        return NULL;
    
    if (NULL == matin)  return NULL;
    if (not_doublematrix(matin)) return NULL;
    if (not_boolmatrix(maskin)) return NULL;
    
    ny = dims[0] = PyArray_DIM(matin, 0);
    nx = dims[1] = PyArray_DIM(matin, 1);
    
    matout = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_DOUBLE);
    
    cin  = PyArray_DATA(matin);
    cout = PyArray_DATA(matout);
    mask = PyArray_DATA(maskin);
    
    buf = malloc(wx*wy*sizeof(*buf));
    offs_x = wx/2; offs_y = wy/2;
    
    for (i = 0; i < ny; i+=step) {
        for (j = 0; j < nx; j+=step) {
	    //if (mask[i*nx + j]) continue;
            for (n = 0, y = i - offs_y; y < i-offs_y+wy; y++) {
                for (x = j - offs_x; x < j-offs_x+wx; x++) {
                    if (x < 0 || x >= nx || y < 0 || y >= ny) continue;
		    idx = y*nx+x;
                    if (mask[idx]) continue;
                    buf[n] = cin[idx]; n++;
                }
            }
            if (n < 3) {
	      cout[i*nx + j] = NAN;
	    } else {
	      cout[i*nx + j] = nr_select(n/2+1, n, buf-1);
	    } 
        }
    }
    
    free(buf);
    return PyArray_Return(matout);
}

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "medfilt",
    NULL,
    -1,
    _medfilt,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_medfilt(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    import_array();
    return m;
}

#else

PyMODINIT_FUNC initmedfilt(void)
{
    PyObject *m;
        
    m = Py_InitModule("medfilt", _medfilt);
    if (m == NULL) {
        return;
    }
    import_array(); 
}
#endif
