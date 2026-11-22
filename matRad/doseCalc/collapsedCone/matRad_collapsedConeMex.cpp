/*
 * matRad_collapsedConeMex: collapsed cone convolution transport kernel
 *
 * call:
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att)
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att, origin)
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att, origin, numThreads)
 *
 * input:
 *   terma:       3D double array of TERMA
 *   rho:         3D double array of relative electron densities (same size)
 *   voxelSize:   voxel size along the three array dimensions [mm]
 *   dirs:        3 x nDirs unit transport directions in array axes
 *   weights:     nComponents x nDirs kernel weights of the exponential
 *                components (e.g. primary and scatter) per direction
 *   att:         nComponents x nDirs attenuation coefficients [1/mm of
 *                water]
 *   origin:      (optional) zero based index of the first voxel of the
 *                arrays in a global grid. The transport lattice is anchored
 *                to the global grid, so that results of sub-volumes of the
 *                same grid are consistent. Default: [0 0 0]
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
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const long maxComponents = 8;

struct Component {
    double weight;
    double att;
};

struct Grid {
    long n[3];
    double s[3];
    long origin[3];
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

inline long roundShift(double x)
{
    return static_cast<long>(std::floor(x + 0.5));
}

void transportDirection(const double* terma, const double* rho, double* dose, const Grid& grid, const double dir[3],
                        const Component* comps, long numComps)
{
    double flow[maxComponents];
    int k = 0;
    double best = -1.0;
    for (int ax = 0; ax < 3; ++ax) {
        const double r = std::fabs(dir[ax]) / grid.s[ax];
        if (r > best) {
            best = r;
            k = ax;
        }
    }
    const int t1 = (k + 1) % 3;
    const int t2 = (k + 2) % 3;
    const long stride[3] = {1, grid.n[0], grid.n[0] * grid.n[1]};
    const long nk = grid.n[k];
    const long n1 = grid.n[t1];
    const long n2 = grid.n[t2];

    const double step = grid.s[k] / std::fabs(dir[k]);
    const double sigma1 = dir[t1] * grid.s[k] / (dir[k] * grid.s[t1]);
    const double sigma2 = dir[t2] * grid.s[k] / (dir[k] * grid.s[t2]);

    std::vector<long> sh1(nk);
    std::vector<long> sh2(nk);
    for (long p = 0; p < nk; ++p) {
        const double globalPlane = static_cast<double>(p + grid.origin[k]);
        sh1[p] = roundShift(globalPlane * sigma1) - grid.origin[t1];
        sh2[p] = roundShift(globalPlane * sigma2) - grid.origin[t2];
    }
    const long minLine1 = -*std::max_element(sh1.begin(), sh1.end());
    const long maxLine1 = n1 - 1 - *std::min_element(sh1.begin(), sh1.end());
    const long minLine2 = -*std::max_element(sh2.begin(), sh2.end());
    const long maxLine2 = n2 - 1 - *std::min_element(sh2.begin(), sh2.end());

    const bool forward = dir[k] > 0.0;

    for (long j0 = minLine2; j0 <= maxLine2; ++j0) {
        for (long i0 = minLine1; i0 <= maxLine1; ++i0) {
            for (long c = 0; c < numComps; ++c) {
                flow[c] = 0.0;
            }
            for (long q = 0; q < nk; ++q) {
                const long p = forward ? q : nk - 1 - q;
                const long i = i0 + sh1[p];
                const long j = j0 + sh2[p];
                if (i < 0 || i >= n1 || j < 0 || j >= n2) {
                    continue;
                }
                const long ix = p * stride[k] + i * stride[t1] + j * stride[t2];
                const double dt = step * rho[ix];
                const double t = terma[ix];
                double d = 0.0;
                for (long c = 0; c < numComps; ++c) {
                    depositStep(dt, t, comps[c], flow[c], d);
                }
                dose[ix] += d;
            }
        }
    }
}

void directionComponents(const double* weights, const double* att, long d, long numComps, Component* comps)
{
    for (long c = 0; c < numComps; ++c) {
        comps[c].weight = weights[d * numComps + c];
        comps[c].att = att[d * numComps + c];
    }
}

void transportAll(const double* terma, const double* rho, double* dose, const Grid& grid, const double* dirs,
                  const double* weights, const double* att, long numComps, long numDirs, std::size_t numVox,
                  int numThreads)
{
#ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = omp_get_max_threads();
    }
    numThreads = static_cast<int>(std::min<long>(numThreads, std::max<long>(numDirs, 1)));
    if (numThreads > 1) {
        std::vector<std::vector<double> > buffers(numThreads);
#pragma omp parallel num_threads(numThreads)
        {
            const int tid = omp_get_thread_num();
            buffers[tid].assign(numVox, 0.0);
            Component comps[maxComponents];
#pragma omp for schedule(static)
            for (long d = 0; d < numDirs; ++d) {
                directionComponents(weights, att, d, numComps, comps);
                transportDirection(terma, rho, buffers[tid].data(), grid, dirs + 3 * d, comps, numComps);
            }
        }
        for (int t = 0; t < numThreads; ++t) {
            if (buffers[t].empty()) {
                continue;
            }
            for (std::size_t v = 0; v < numVox; ++v) {
                dose[v] += buffers[t][v];
            }
        }
        return;
    }
#else
    (void)numThreads;
#endif
    Component comps[maxComponents];
    for (long d = 0; d < numDirs; ++d) {
        directionComponents(weights, att, d, numComps, comps);
        transportDirection(terma, rho, dose, grid, dirs + 3 * d, comps, numComps);
    }
}

}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs < 6 || nrhs > 8) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:nrhs", "Six to eight inputs required.");
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
    Grid grid = {{1, 1, 1}, {0.0, 0.0, 0.0}, {0, 0, 0}};
    for (mwSize d = 0; d < numDims; ++d) {
        grid.n[d] = static_cast<long>(dims[d]);
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
    const long numDirs = static_cast<long>(mxGetN(prhs[3]));
    const long numComps = static_cast<long>(mxGetM(prhs[4]));
    if (numComps < 1 || numComps > maxComponents) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:components", "Between one and eight kernel components are supported.");
    }
    if (static_cast<long>(mxGetN(prhs[4])) != numDirs || mxGetM(prhs[5]) != mxGetM(prhs[4]) ||
        mxGetN(prhs[5]) != mxGetN(prhs[4])) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:kernel",
                          "Kernel weights and attenuation coefficients must be nComponents x nDirs arrays.");
    }

    const double* s = mxGetPr(prhs[2]);
    for (int d = 0; d < 3; ++d) {
        grid.s[d] = s[d];
    }
    if (nrhs > 6 && !mxIsEmpty(prhs[6])) {
        if (mxGetNumberOfElements(prhs[6]) != 3) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:origin", "Origin must have three elements.");
        }
        const double* o = mxGetPr(prhs[6]);
        for (int d = 0; d < 3; ++d) {
            grid.origin[d] = static_cast<long>(o[d]);
        }
    }
    int numThreads = 0;
    if (nrhs > 7) {
        numThreads = static_cast<int>(mxGetScalar(prhs[7]));
    }

    plhs[0] = mxCreateNumericArray(numDims, dims, mxDOUBLE_CLASS, mxREAL);
    transportAll(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(plhs[0]), grid, mxGetPr(prhs[3]), mxGetPr(prhs[4]),
                 mxGetPr(prhs[5]), numComps, numDirs, numVox, numThreads);
}
