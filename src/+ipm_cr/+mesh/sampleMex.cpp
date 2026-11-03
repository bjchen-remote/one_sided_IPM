#include "mex.h"
#include <cmath>

// CSC traversal of a sparse query map. Each output starts at its strongest
// local donor and accumulates weighted differences. No transpose, expanded
// per-nonzero difference arrays or dense query-by-donor storage are needed.
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3 || nlhs != 1 || !mxIsSparse(prhs[0]) ||
        !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !mxIsDouble(prhs[1]) || mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) ||
        !mxIsDouble(prhs[2]) || mxIsSparse(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgIdAndTxt("ipm_cr:SampleNative", "Use a real sparse map, double donors and anchor indices.");
    const mwSize m = mxGetM(prhs[0]), n = mxGetN(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != n || mxGetNumberOfElements(prhs[2]) != m)
        mexErrMsgIdAndTxt("ipm_cr:SampleNative", "Map, donor and anchor dimensions differ.");
    const double* f = mxGetDoubles(prhs[1]);
    const double* a = mxGetDoubles(prhs[2]);
    const double* w = mxGetDoubles(prhs[0]);
    const mwIndex* row = mxGetIr(prhs[0]);
    const mwIndex* start = mxGetJc(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double* value = mxGetDoubles(plhs[0]);
    for (mwIndex i = 0; i < m; ++i) {
        if (!std::isfinite(a[i]) || a[i] < 0 || a[i] > n || std::floor(a[i]) != a[i])
            mexErrMsgIdAndTxt("ipm_cr:SampleNative", "An anchor index is invalid.");
    }
    bool constant = n > 0;
    for (mwIndex j = 1; j < n && constant; ++j) constant = f[j] == f[0];
    if (constant) {
        for (mwIndex i = 0; i < m; ++i) if (a[i] > 0) value[i] = f[0];
        return;
    }
    for (mwIndex j = 0; j < n; ++j)
        for (mwIndex k = start[j]; k < start[j + 1]; ++k) {
            const mwIndex i = row[k];
            if (a[i] <= 0)
                mexErrMsgIdAndTxt("ipm_cr:SampleNative", "A nonempty query has no local anchor.");
            value[i] += w[k] * (f[j] - f[static_cast<mwIndex>(a[i]) - 1]);
        }
    for (mwIndex i = 0; i < m; ++i)
        if (a[i] > 0) value[i] += f[static_cast<mwIndex>(a[i]) - 1];
}
