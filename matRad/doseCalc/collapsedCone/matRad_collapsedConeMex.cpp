/*
 * matRad_collapsedConeMex: collapsed cone convolution transport kernel
 *
 * call:
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, A, a, B, b)
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, A, a, B, b, numThreads)
 *
 * input:
 *   terma:       3D double array of TERMA
 *   rho:         3D double array of relative electron densities (same size)
 *   voxelSize:   voxel size along the three array dimensions [mm]
 *   dirs:        3 x nDirs unit transport directions in array axes
 *   A, a, B, b:  1 x nDirs kernel weights and attenuation coefficients
 *                [1/mm of water] per direction
 *   numThreads:  (optional) number of OpenMP threads, 0 = default
 *
 * output:
 *   dose:        3D double array of deposited dose (same size as terma)
 *
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 * Copyright 2026 the matRad development team.
 *
 * This file is part of the matRad project. It is subject to the license
 * terms in the LICENSE file found in the top-level directory of this
 * distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
 * of the matRad project, including this file, may be copied, modified,
 * propagated, or distributed except according to the terms contained in the
 * LICENSE file.
 *
 * %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 */

#include "mex.h"
#include <cmath>
#include <cstddef>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Component {
    double weight;
    double att;
};

inline void depositStep(double dt, double t, const Component& c, double& p, double& d)
{
    if (c.weight == 0.0) {
        return;
    }
    const double x = c.att * dt;
    double g;
    double oneMinusE;
    if (x > 0.0) {
        oneMinusE = -std::expm1(-x);
        g = oneMinusE / x;
    } else {
        oneMinusE = 0.0;
        g = 1.0;
    }
    const double source = t * c.weight;
    d += p * g + source * (1.0 - g);
    p = p * (1.0 - oneMinusE) + source * oneMinusE;
}

inline long shiftAt(long q, double delta)
{
    return static_cast<long>(std::floor(static_cast<double>(q) * delta + 0.5));
}

void transportDirection(const double* terma, const double* rho, double* dose, const long n[3], const double s[3],
                        const double dir[3], const Component& ca, const Component& cb)
{
    int k = 0;
    double best = -1.0;
    for (int ax = 0; ax < 3; ++ax) {
        const double r = std::fabs(dir[ax]) / s[ax];
        if (r > best) {
            best = r;
            k = ax;
        }
    }
    const int t1 = (k + 1) % 3;
    const int t2 = (k + 2) % 3;
    const long stride[3] = {1, n[0], n[0] * n[1]};

    const double step = s[k] / std::fabs(dir[k]);
    const double delta1 = dir[t1] * s[k] / (std::fabs(dir[k]) * s[t1]);
    const double delta2 = dir[t2] * s[k] / (std::fabs(dir[k]) * s[t2]);
    const bool forward = dir[k] > 0.0;
    const long nk = n[k];

    const long last1 = shiftAt(nk - 1, delta1);
    const long last2 = shiftAt(nk - 1, delta2);
    const long min1 = -(last1 > 0 ? last1 : 0);
    const long max1 = n[t1] - 1 - (last1 < 0 ? last1 : 0);
    const long min2 = -(last2 > 0 ? last2 : 0);
    const long max2 = n[t2] - 1 - (last2 < 0 ? last2 : 0);
    const long numLines1 = max1 - min1 + 1;
    const long numLines2 = max2 - min2 + 1;
    const long numLines = numLines1 * numLines2;

    std::vector<long> sh1(nk);
    std::vector<long> sh2(nk);
    for (long q = 0; q < nk; ++q) {
        sh1[q] = shiftAt(q, delta1);
        sh2[q] = shiftAt(q, delta2);
    }

#pragma omp parallel for schedule(static)
    for (long line = 0; line < numLines; ++line) {
        const long i0 = min1 + line % numLines1;
        const long j0 = min2 + line / numLines1;
        double pa = 0.0;
        double pb = 0.0;
        for (long q = 0; q < nk; ++q) {
            const long i = i0 + sh1[q];
            const long j = j0 + sh2[q];
            if (i < 0 || i >= n[t1] || j < 0 || j >= n[t2]) {
                continue;
            }
            const long plane = forward ? q : nk - 1 - q;
            const long ix = plane * stride[k] + i * stride[t1] + j * stride[t2];
            const double dt = step * rho[ix];
            const double t = terma[ix];
            double d = 0.0;
            depositStep(dt, t, ca, pa, d);
            depositStep(dt, t, cb, pb, d);
            dose[ix] += d;
        }
    }
}

}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 8 || nrhs > 9) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:nrhs", "Eight or nine inputs required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:nlhs", "One output allowed.");
    }
    for (int i = 0; i < nrhs; ++i) {
        if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i])) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:type", "All inputs must be real, full double arrays.");
        }
    }

    const mwSize numDims = mxGetNumberOfDimensions(prhs[0]);
    const mwSize* dims = mxGetDimensions(prhs[0]);
    if (numDims > 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:dims", "TERMA must be a 3D array.");
    }
    long n[3] = {1, 1, 1};
    for (mwSize d = 0; d < numDims; ++d) {
        n[d] = static_cast<long>(dims[d]);
    }
    const std::size_t numVox = mxGetNumberOfElements(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != numVox) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:size", "TERMA and density must have the same size.");
    }
    if (mxGetNumberOfElements(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:voxelSize", "Voxel size must have three elements.");
    }
    if (mxGetM(prhs[3]) != 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:dirs", "Directions must be a 3 x nDirs array.");
    }
    const std::size_t numDirs = mxGetN(prhs[3]);
    for (int i = 4; i < 8; ++i) {
        if (mxGetNumberOfElements(prhs[i]) != numDirs) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:kernel", "Kernel parameters must have one value per direction.");
        }
    }

    const double* terma = mxGetPr(prhs[0]);
    const double* rho = mxGetPr(prhs[1]);
    const double* s = mxGetPr(prhs[2]);
    const double* dirs = mxGetPr(prhs[3]);
    const double* wA = mxGetPr(prhs[4]);
    const double* aA = mxGetPr(prhs[5]);
    const double* wB = mxGetPr(prhs[6]);
    const double* aB = mxGetPr(prhs[7]);

#ifdef _OPENMP
    if (nrhs > 8) {
        const int numThreads = static_cast<int>(mxGetScalar(prhs[8]));
        if (numThreads > 0) {
            omp_set_num_threads(numThreads);
        }
    }
#endif

    plhs[0] = mxCreateNumericArray(numDims, dims, mxDOUBLE_CLASS, mxREAL);
    double* dose = mxGetPr(plhs[0]);
    const double voxelSize[3] = {s[0], s[1], s[2]};

    for (std::size_t d = 0; d < numDirs; ++d) {
        const Component ca = {wA[d], aA[d]};
        const Component cb = {wB[d], aB[d]};
        transportDirection(terma, rho, dose, n, voxelSize, dirs + 3 * d, ca, cb);
    }
}
