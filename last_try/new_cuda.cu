#ifndef _MCM_CTLDEV_DEBUGCORR_H
#define _MCM_CTLDEV_DEBUGCORR_H

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

struct ampl_and_phase_debug {
    double ampl;
    double phase;    // In radians.
};

__global__ void phase_rotation_kernel(const double* decodeddata, double curfreq, double curphase, double* re, double* im, int datablockSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < datablockSize) {
        double decodedval = decodeddata[idx];
        double phi_i = curphase + curfreq * (double)idx;
        phi_i = (phi_i - floor(phi_i)) * 2.0 * M_PI;
        re[idx] = decodedval * cos(phi_i);
        im[idx] = decodedval * sin(phi_i);
    }
}

__global__ void low_pass_filter_kernel(const double* re, const double* im, double* reLPF, double* imLPF, int const_LPFThreshold, int datablockSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < datablockSize) {
        reLPF[idx] = imLPF[idx] = 0.0;

        int lastDataunitIndex = idx + const_LPFThreshold;
        if (lastDataunitIndex > 8192)
            lastDataunitIndex = 8192;

        for (int k = idx; k < lastDataunitIndex; k++) {
            reLPF[idx] += re[k];
            imLPF[idx] += im[k];
        }

        reLPF[idx] /= (double)const_LPFThreshold;
        imLPF[idx] /= (double)const_LPFThreshold;
    }
}

__global__ void compute_sum_kernel(const double* reLPF, const double* imLPF, double* reSum, double* imSum, int datablockSize) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < datablockSize) {
        double re2 = reLPF[idx] * reLPF[idx] - imLPF[idx] * imLPF[idx];
        double im2 = 2.0 * reLPF[idx] * imLPF[idx];

        double re4 = re2 * re2 - im2 * im2;
        double im4 = 2.0 * re2 * im2;

        atomicAdd(reSum, re4);
        atomicAdd(imSum, im4);
    }
}

void controller_debugcorr(const double* origdata, double curfreq, double curphase, struct ampl_and_phase_debug* outdata) {
    const int datablockSize = 8192;
    const int threadsPerBlock = 256;
    const int numBlocks = (datablockSize + threadsPerBlock - 1) / threadsPerBlock;

    double* decodeddata;
    double* re;
    double* im;
    double* reLPF;
    double* imLPF;
    double* reSum;
    double* imSum;

    cudaMalloc((void**)&decodeddata, datablockSize * sizeof(double));
    cudaMalloc((void**)&re, datablockSize * sizeof(double));
    cudaMalloc((void**)&im, datablockSize * sizeof(double));
    cudaMalloc((void**)&reLPF, datablockSize * sizeof(double));
    cudaMalloc((void**)&imLPF, datablockSize * sizeof(double));
    cudaMalloc((void**)&reSum, sizeof(double));
    cudaMalloc((void**)&imSum, sizeof(double));

    cudaMemcpy(decodeddata, origdata, datablockSize * sizeof(double), cudaMemcpyHostToDevice);

    phase_rotation_kernel<<<numBlocks, threadsPerBlock>>>(decodeddata, curfreq, curphase, re, im, datablockSize);
    cudaDeviceSynchronize();

    low_pass_filter_kernel<<<numBlocks, threadsPerBlock>>>(re, im, reLPF, imLPF, 10, datablockSize);
    cudaDeviceSynchronize();

    compute_sum_kernel<<<numBlocks, threadsPerBlock>>>(reLPF, imLPF, reSum, imSum, datablockSize);
    cudaDeviceSynchronize();

    double host_reSum, host_imSum;
    cudaMemcpy(&host_reSum, reSum, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_imSum, imSum, sizeof(double), cudaMemcpyDeviceToHost);

    host_reSum *= 1000.0 / datablockSize;
    host_imSum *= 1000.0 / datablockSize;

    outdata->ampl = sqrt(host_reSum * host_reSum + host_imSum * host_imSum);
    outdata->phase = atan2(host_imSum, host_reSum);

    cudaFree(decodeddata);
    cudaFree(re);
    cudaFree(im);
    cudaFree(reLPF);
    cudaFree(imLPF);
    cudaFree(reSum);
    cudaFree(imSum);
}

#endif /* of #ifndef _MCM_CTLDEV_DEBUGCORR_H */
