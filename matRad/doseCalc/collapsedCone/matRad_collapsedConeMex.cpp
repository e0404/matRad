/*
 * matRad_collapsedConeMex: collapsed cone convolution transport kernel
 *
 * call:
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att)
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att, origin)
 *   dose = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, weights, att, origin, numThreads)
 *
 * input:
 *   terma:       3D single or double array of TERMA
 *   rho:         3D single or double array of relative electron densities (same size and class as terma)
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
 *   dose:        3D array of deposited dose (same size and class as terma)
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
#include <cstring>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const std::ptrdiff_t maxComponents = 8;

template <typename T>
struct Component {
    T weight;
    T att;
};

template <typename T>
struct Grid {
    std::ptrdiff_t n[3];
    T s[3];
    std::ptrdiff_t origin[3];
};

template <typename T>
inline void depositStep(T dt, T t, const Component<T>& c, T& p, T& d)
{
    if (c.weight == T(0)) {
        return;
    }
    const T x = c.att * dt;
    T g;
    T oneMinusE;
    if (x > T(0)) {
        oneMinusE = -std::expm1(-x);
        g = oneMinusE / x;
    } else {
        oneMinusE = T(0);
        g = T(1);
    }
    const T source = t * c.weight;
    d += p * g + source * (T(1) - g);
    p = p * (T(1) - oneMinusE) + source * oneMinusE;
}

template <typename T>
inline std::ptrdiff_t roundShift(T x)
{
    return static_cast<std::ptrdiff_t>(std::floor(x + T(0.5)));
}

template <typename T>
void transportDirection(const T* terma, const T* rho, T* dose, const Grid<T>& grid, const T dir[3],
                        const Component<T>* comps, std::ptrdiff_t numComps)
{
    T flow[maxComponents];
    int k = 0;
    T best = T(-1);
    for (int ax = 0; ax < 3; ++ax) {
        const T r = std::fabs(dir[ax]) / grid.s[ax];
        if (r > best) {
            best = r;
            k = ax;
        }
    }
    const int t1 = (k + 1) % 3;
    const int t2 = (k + 2) % 3;
    const std::ptrdiff_t stride[3] = {1, grid.n[0], grid.n[0] * grid.n[1]};
    const std::ptrdiff_t nk = grid.n[k];
    const std::ptrdiff_t n1 = grid.n[t1];
    const std::ptrdiff_t n2 = grid.n[t2];

    const T step = grid.s[k] / std::fabs(dir[k]);
    const T sigma1 = dir[t1] * grid.s[k] / (dir[k] * grid.s[t1]);
    const T sigma2 = dir[t2] * grid.s[k] / (dir[k] * grid.s[t2]);

    std::vector<std::ptrdiff_t> sh1(nk);
    std::vector<std::ptrdiff_t> sh2(nk);
    for (std::ptrdiff_t p = 0; p < nk; ++p) {
        const T globalPlane = static_cast<T>(p + grid.origin[k]);
        sh1[p] = roundShift(globalPlane * sigma1) - grid.origin[t1];
        sh2[p] = roundShift(globalPlane * sigma2) - grid.origin[t2];
    }
    const std::ptrdiff_t minLine1 = -*std::max_element(sh1.begin(), sh1.end());
    const std::ptrdiff_t maxLine1 = n1 - 1 - *std::min_element(sh1.begin(), sh1.end());
    const std::ptrdiff_t minLine2 = -*std::max_element(sh2.begin(), sh2.end());
    const std::ptrdiff_t maxLine2 = n2 - 1 - *std::min_element(sh2.begin(), sh2.end());

    const bool forward = dir[k] > T(0);

    for (std::ptrdiff_t j0 = minLine2; j0 <= maxLine2; ++j0) {
        for (std::ptrdiff_t i0 = minLine1; i0 <= maxLine1; ++i0) {
            for (std::ptrdiff_t c = 0; c < numComps; ++c) {
                flow[c] = T(0);
            }
            for (std::ptrdiff_t q = 0; q < nk; ++q) {
                const std::ptrdiff_t p = forward ? q : nk - 1 - q;
                const std::ptrdiff_t i = i0 + sh1[p];
                const std::ptrdiff_t j = j0 + sh2[p];
                if (i < 0 || i >= n1 || j < 0 || j >= n2) {
                    continue;
                }
                const std::ptrdiff_t ix = p * stride[k] + i * stride[t1] + j * stride[t2];
                const T t = terma[ix];
                if (t == T(0)) {
                    bool zeroFlow = true;
                    for (std::ptrdiff_t c = 0; c < numComps; ++c) {
                        if (flow[c] != T(0)) {
                            zeroFlow = false;
                            break;
                        }
                    }
                    if (zeroFlow) {
                        continue;
                    }
                }
                const T dt = step * rho[ix];
                T d = T(0);
                for (std::ptrdiff_t c = 0; c < numComps; ++c) {
                    depositStep(dt, t, comps[c], flow[c], d);
                }
                dose[ix] += d;
            }
        }
    }
}

template <typename T>
void directionComponents(const T* weights, const T* att, std::ptrdiff_t d, std::ptrdiff_t numComps,
                         Component<T>* comps)
{
    for (std::ptrdiff_t c = 0; c < numComps; ++c) {
        comps[c].weight = weights[d * numComps + c];
        comps[c].att = att[d * numComps + c];
    }
}

template <typename T>
void transportAll(const T* terma, const T* rho, T* dose, const Grid<T>& grid, const T* dirs,
                  const T* weights, const T* att, std::ptrdiff_t numComps, std::ptrdiff_t numDirs,
                  std::size_t numVox, int numThreads)
{
#ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = omp_get_max_threads();
    }
    numThreads = static_cast<int>(std::min<std::ptrdiff_t>(numThreads, std::max<std::ptrdiff_t>(numDirs, 1)));
    if (numThreads > 1) {
        std::vector<std::vector<T> > buffers(numThreads);
#pragma omp parallel num_threads(numThreads)
        {
            const int tid = omp_get_thread_num();
            buffers[tid].assign(numVox, T(0));
            Component<T> comps[maxComponents];
#pragma omp for schedule(static)
            for (std::ptrdiff_t d = 0; d < numDirs; ++d) {
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
    Component<T> comps[maxComponents];
    for (std::ptrdiff_t d = 0; d < numDirs; ++d) {
        directionComponents(weights, att, d, numComps, comps);
        transportDirection(terma, rho, dose, grid, dirs + 3 * d, comps, numComps);
    }
}

template <typename T>
std::vector<T> readNumericVector(const mxArray* arr)
{
    const std::size_t n = mxGetNumberOfElements(arr);
    std::vector<T> out(n);
    if (mxIsDouble(arr)) {
        const double* ptr = static_cast<const double*>(mxGetData(arr));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = static_cast<T>(ptr[i]);
        }
    } else {
        const float* ptr = static_cast<const float*>(mxGetData(arr));
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = static_cast<T>(ptr[i]);
        }
    }
    return out;
}

template <typename T>
void runMex(int nrhs, const mxArray* prhs[], mxArray* plhs[], mwSize numDims, const mwSize* dims,
            std::size_t numVox, std::ptrdiff_t numDirs, std::ptrdiff_t numComps, int numThreads)
{
    Grid<T> grid = {{1, 1, 1}, {T(0), T(0), T(0)}, {0, 0, 0}};
    for (mwSize d = 0; d < numDims; ++d) {
        grid.n[d] = static_cast<std::ptrdiff_t>(dims[d]);
    }
    const std::vector<T> s = readNumericVector<T>(prhs[2]);
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(s[d]) || s[d] <= T(0)) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:voxelSize",
                              "Voxel size elements must be finite and strictly positive.");
        }
        grid.s[d] = s[d];
    }
    if (nrhs > 6 && !mxIsEmpty(prhs[6])) {
        const std::vector<double> o = readNumericVector<double>(prhs[6]);
        for (int d = 0; d < 3; ++d) {
            grid.origin[d] = static_cast<std::ptrdiff_t>(o[d]);
        }
    }

    const std::vector<T> dirs = readNumericVector<T>(prhs[3]);
    for (std::ptrdiff_t d = 0; d < numDirs; ++d) {
        const T dx = dirs[3 * d];
        const T dy = dirs[3 * d + 1];
        const T dz = dirs[3 * d + 2];
        const T norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (!std::isfinite(norm) || std::fabs(norm - T(1)) > T(1e-3)) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:dirs",
                              "Each direction vector must be finite and unit-length.");
        }
    }

    const std::vector<T> weights = readNumericVector<T>(prhs[4]);
    const std::vector<T> att = readNumericVector<T>(prhs[5]);

    const T* termaPtr = static_cast<const T*>(mxGetData(prhs[0]));
    const T* rhoPtr = static_cast<const T*>(mxGetData(prhs[1]));
    T* dosePtr = static_cast<T*>(mxGetData(plhs[0]));

    transportAll<T>(termaPtr, rhoPtr, dosePtr, grid, dirs.data(), weights.data(), att.data(),
                    numComps, numDirs, numVox, numThreads);
}

}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs == 1 && mxIsChar(prhs[0])) {
        char cmd[16];
        if (mxGetString(prhs[0], cmd, sizeof(cmd)) == 0 && std::strcmp(cmd, "unlock") == 0) {
            if (mexIsLocked()) {
                mexUnlock();
            }
            return;
        }
    }

#if defined(_MSC_VER) && defined(_OPENMP)
    if (!mexIsLocked()) {
        mexLock();
    }
#endif

    if (nrhs < 6 || nrhs > 8) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:nrhs", "Six to eight inputs required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:nlhs", "One output allowed.");
    }
    for (int i = 0; i < nrhs; ++i) {
        if ((!mxIsDouble(prhs[i]) && !mxIsSingle(prhs[i])) || mxIsComplex(prhs[i]) || mxIsSparse(prhs[i])) {
            mexErrMsgIdAndTxt("matRad:collapsedConeMex:type", "All inputs must be real, full single or double arrays.");
        }
    }

    const bool isSingle = mxIsSingle(prhs[0]);
    if (mxIsSingle(prhs[1]) != isSingle) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:precisionMismatch", "TERMA and density must have the same precision.");
    }

    const mwSize numDims = mxGetNumberOfDimensions(prhs[0]);
    const mwSize* dims = mxGetDimensions(prhs[0]);
    if (numDims > 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:dims", "TERMA must be a 3D array.");
    }
    const std::size_t numVox = mxGetNumberOfElements(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != numVox) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:size", "TERMA and density must have the same size.");
    }
    if (mxGetNumberOfElements(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:voxelSize", "Voxel size must have three elements.");
    }
    if (mxGetM(prhs[3]) != 3 || mxGetN(prhs[3]) < 1) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:dirs", "Directions must be a 3 x nDirs array with at least one direction.");
    }
    const std::ptrdiff_t numDirs = static_cast<std::ptrdiff_t>(mxGetN(prhs[3]));
    const std::ptrdiff_t numComps = static_cast<std::ptrdiff_t>(mxGetM(prhs[4]));
    if (numComps < 1 || numComps > maxComponents) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:components", "Between one and eight kernel components are supported.");
    }
    if (static_cast<std::ptrdiff_t>(mxGetN(prhs[4])) != numDirs || mxGetM(prhs[5]) != mxGetM(prhs[4]) ||
        mxGetN(prhs[5]) != mxGetN(prhs[4])) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:kernel",
                          "Kernel weights and attenuation coefficients must be nComponents x nDirs arrays.");
    }
    if (nrhs > 6 && !mxIsEmpty(prhs[6]) && mxGetNumberOfElements(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("matRad:collapsedConeMex:origin", "Origin must have three elements.");
    }

    int numThreads = 0;
    if (nrhs > 7 && !mxIsEmpty(prhs[7])) {
        numThreads = static_cast<int>(mxGetScalar(prhs[7]));
    }

    if (isSingle) {
        plhs[0] = mxCreateNumericArray(numDims, dims, mxSINGLE_CLASS, mxREAL);
        runMex<float>(nrhs, prhs, plhs, numDims, dims, numVox, numDirs, numComps, numThreads);
    } else {
        plhs[0] = mxCreateNumericArray(numDims, dims, mxDOUBLE_CLASS, mxREAL);
        runMex<double>(nrhs, prhs, plhs, numDims, dims, numVox, numDirs, numComps, numThreads);
    }
}
