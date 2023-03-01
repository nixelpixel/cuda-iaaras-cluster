#include <stdio.h>
#include "iostream"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "chrono"

//#include "code.cu"

#if defined(_WIN32) || defined(_WIN64)
#	include "fftw3.h"
#else

#	include <fftw3.h>

#endif

#define COUNT 4

#define THREADS_GPU_COUNT 5
#define ARR_SIZE 1000000

#define POW2 8192
//#define POW2 10000

struct DFHeader {
    char hat[16];
};

struct d{
    double r;
    double i;
};

struct Frequency_pair {

    double C;
    double S;
};

struct DF {
    struct DFHeader header;
    char data[10000];
    float decoded_data[10000];
    struct Frequency_pair phase_data[10000];
    struct Frequency_pair after_filter_data[10000];
    struct Frequency_pair fourfold_phase_data[10000];
    double Cs, Ss;
    double ampl;
    double phi;
    double f_next;


};

struct Decoder {

    FILE *file;
    struct DF * df;

};

__global__ void add(int a, int b, int *c)
{
    *c = a + b;
}

void getInfo(){
    int count;
    cudaGetDeviceCount(&count);
    cudaDeviceProp prop{};
    for (int i = 0; i<count; i++){
        cudaGetDeviceProperties(&prop, i);
        printf("--- Общая информация об устройстве %d ---\n", i);
        printf("Имя: %s\n", prop.name);
        printf("Вычислительные возможности: %d.%d\n", prop.major, prop.minor);
        printf("Тактовая частота: %d\n", prop.clockRate);
        printf("Перекрытие копирования: ");
        if(prop.deviceOverlap){
            printf("Разрешено\n");
        } else {
            printf("Запрещено");
        }
        printf("Тайм-аут выполнения ядра: ");
        if (prop.kernelExecTimeoutEnabled){
            printf("Включен\n");
        } else {
            printf("Выключен\n");
        }
        printf("--- Информация о памяти для устройства %d ---\n", i);
        printf("Всего глобальной памяти: %ld\n", prop.totalGlobalMem);
        printf("Всего константной памяти: %ld\n", prop.totalConstMem);
        printf("Максимальный шаг: %ld\n", prop.memPitch);
        printf("Выравнивание текстур: %ld\n", prop.textureAlignment);

        printf("--- Информация о мультипроцессорах для устройства %d ---\n", i);
        printf("Количество мультипроцессоров: %d\n", prop.multiProcessorCount);
        printf("Разделяемая память на один МП: %ld\n", prop.sharedMemPerBlock);
        printf("Регистров на один МП: %d\n", prop.regsPerBlock);
        printf("Нитей в варпе: %d\n", prop.warpSize);
        printf("Максимальное количество нитей в блоке: %d\n", prop.maxThreadsPerBlock);
        printf("Макс. количество нитей по измерениям: (%d, %d, %d)\n",prop.maxThreadsDim[0],prop.maxThreadsDim[1],prop.maxThreadsDim[2]);
        printf("Максимальные размеры сетки: (%d, %d, %d)\n",prop.maxGridSize[0],prop.maxGridSize[1],prop.maxGridSize[2]);
        printf("\n");
    }
}

void decoding(char in, float * res){

    char arr;

    arr = (in >> 2) & 0x3;

    if (arr == 0x0){
        *res = 0.3333;
    }
    if (arr == 0x1){
        *res = -1.0;
    }
    if (arr == 0x2){
        *res = -0.3333;
    }
    if (arr == 0x3){
        *res = -1.0;
    }

}

void read_file(char *file_name, struct Decoder *decoder) {

    decoder->file = fopen(file_name, "rb");
    if (!decoder->file) {
        perror("File error");
        return;
    }

    FILE *f = decoder->file;
    decoder->df = (struct DF *) malloc(sizeof(struct DF) * COUNT);
    for (int i = 0; i < COUNT; ++i) {
        fread(decoder->df[i].header.hat, 16, 1, f);
        fread(decoder->df[i].data, 10000, 1, f);
//
//        for (int j = 0; j < 16; ++j) {
////            printf("[%d] = %x\n", j, decoder->df[i].header.hat[j]);
////            printf("[%d] data = %x\n", j, decoder->df[i].data[j]);
//        }
//        printf("\n");
    }

    //DECODING
    for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < 10000; ++j) {
            decoding(decoder->df[i].data[j],&decoder->df[i].decoded_data[j]);
//            printf("data[%d] = %f\n",i, decoder->df[j].decoded_data[i]);
        }
    }
    fclose(decoder->file);


}


void record_float_to_file(float array[10000] ) {

    FILE *fp;
    fp = fopen("../floats.txt", "w");

    for (unsigned i = 0; i < 10000; i++) {
        fprintf(fp, "%f\n", array[i]);

    }

    fclose(fp);

}

float calc_alpha(long int t, float f){

    double alpha = f * (double)t;
//    phi += M_PI_4;
    alpha-=floor(alpha);
    alpha*=2*M_PI;
    return alpha;
}

struct Frequency_pair calc_freq_pair(long int t, /*double phi,*/ float I){
    struct Frequency_pair res;
    float f = 1851.0 / 8192.0;
    float alpha=calc_alpha(t,f);

    res.C = I * cos(alpha);
    res.S = I * sin(alpha);
    return res;
}

void phase_rotation(struct DF * df) {
    FILE *fp;
    fp = fopen("../phase.txt", "w");
    for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < 10000; ++j) {
            df[i].phase_data[j].S = calc_freq_pair(i * 10000 + j,df[i].decoded_data[j]).S;
            df[i].phase_data[j].C = calc_freq_pair(i * 10000 + j,df[i].decoded_data[j]).C;
//            printf("phase[%d][%d]: C = %f| S = %f\n", i,j,df[i].phase_data[j].C,df[i].phase_data[j].S);

            fprintf(fp, "%f %f\n",df[i].phase_data[j].C,df[i].phase_data[j].S);

        }

    }
    fclose(fp);

}

__global__ void FilterLowerFreqGPU(struct DF * device_df){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

    for(int k = 0; k < COUNT; ++k){
        for (int i = 0; i < 10000; ++i) {
            local_avg_re = 0;
            local_avg_im = 0;

            for (int j = -4; j <= 3; ++j) {
                index = i+j;
                if (index < 0){
                    index = 0;
                }
                else if ( index >= 10000) {
                    index = 9999;
                }
                local_avg_re += device_df[k].phase_data[index].C;
                local_avg_im = device_df[k].phase_data[index].S;
            }

            device_df[k].after_filter_data[i].C = local_avg_re / 8.0;
            device_df[k].after_filter_data[i].S = local_avg_im / 8.0;
        }
    }


}

void FilterLowerFreq(struct DF * df){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

    for(int k = 0; k < COUNT; ++k){
        for (int i = 0; i < 10000; ++i) {
            local_avg_re = 0;
            local_avg_im = 0;

            for (int j = -4; j <= 3; ++j) {
                index = i+j;
                if (index < 0){
                    index = 0;
                }
                else if ( index >= 10000) {
                    index = 9999;
                }
                local_avg_re += df[k].phase_data[index].C;
                local_avg_im = df[k].phase_data[index].S;
            }

            df[k].after_filter_data[i].C = local_avg_re / 8.0;
            df[k].after_filter_data[i].S = local_avg_im / 8.0;
        }
    }


    //фУРЬЕ =============================================

    double after_arr[10000];
    struct d fourier[10000];

    fftw_plan plan;
    plan = fftw_plan_dft_1d(POW2, (fftw_complex *) df->after_filter_data, (fftw_complex *) fourier, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    FILE *file2;
    file2 = fopen("../filter_fourier.txt", "w");

    for (int i = 0; i < POW2; ++i) {
        after_arr[i] = sqrt(fourier[i].r * fourier[i].r + fourier[i].i * fourier[i].i);

        fprintf(file2, "%lf\n", after_arr[i]);
//        printf( "%lf\n", arr[i]);
    }
    fftw_destroy_plan(plan);
    fclose(file2);
//    ====================================================

}



struct Frequency_pair frequency_x4(struct Frequency_pair incoming) {

    struct Frequency_pair frequency_2;
    frequency_2.C = (incoming.C * incoming.C) - (incoming.S * incoming.S);
    frequency_2.S = 2 * incoming.C * incoming.S;

    struct Frequency_pair frequency_4;
    frequency_4.C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
    frequency_4.S = 2 * frequency_2.C * frequency_2.S;

//    return frequency_2;
    return frequency_4;
}

__global__ void frequency_x4_gpu(struct Frequency_pair incoming, struct Frequency_pair out){
    printf("here\n");
    struct Frequency_pair frequency_2;
    frequency_2.C = (incoming.C * incoming.C) - (incoming.S * incoming.S);
    frequency_2.S = 2 * incoming.C * incoming.S;

    out.C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
    out.S = 2 * frequency_2.C * frequency_2.S;
}

void frequency_multiply(struct DF * device_df){
//        for (int i = 0; i < COUNT; ++i) {

            for (int j = 0; j < 10000; ++j) {

                frequency_x4_gpu<<<1,1>>>(device_df[0].after_filter_data[j],device_df[0].fourfold_phase_data[j]);
//                df[i].fourfold_phase_data[j].S = frequency_x4_gpu(df[i].after_filter_data[j]).S;
//                df[i].fourfold_phase_data[j].C = frequency_x4_gpu(df[i].after_filter_data[j]).C;
//                df[i].fourfold_phase_data[j].S = frequency_x4_gpu(df[i].after_filter_data[j]).S;
            }
//    }
}

__global__ void frequency_fourfold_GPU(struct DF * device_df){
    for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < 10000; ++j) {
            struct Frequency_pair frequency_2;
            frequency_2.C = (device_df[i].after_filter_data[j].C * device_df[i].after_filter_data[j].C) - (device_df[i].after_filter_data[j].S * device_df[i].after_filter_data[j].S);
            frequency_2.S = 2 * device_df[i].after_filter_data[j].C * device_df[i].after_filter_data[j].S;

            struct Frequency_pair frequency_4;
            frequency_4.C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
            frequency_4.S = 2 * frequency_2.C * frequency_2.S;

            device_df[i].fourfold_phase_data[j].C = frequency_4.C;
            device_df[i].fourfold_phase_data[j].S = frequency_4.S;
        }
    }
}

void frequency_fourfold(struct DF * df){
//    struct DF * device_df;
//
//    int error_1;
//    error_1 = cudaMalloc((void**) &device_df, sizeof(struct DF) * COUNT);
//    if (error_1){
//        std::cout<<"ERROR_9";
//    }

    for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < 10000; ++j) {
            df[i].fourfold_phase_data[j].C = frequency_x4(df[i].after_filter_data[j]).C;
            df[i].fourfold_phase_data[j].S = frequency_x4(df[i].after_filter_data[j]).S;
        }
    }
//    error_1 = cudaMemcpy(device_df, df, sizeof(struct DF) * COUNT,  cudaMemcpyHostToDevice);
//    if (error_1){
//        std::cout<<"ERROR_10";
//    }
//    frequency_multiply(device_df);
//    error_1 = cudaMemcpy(df, device_df, sizeof(struct DF) * COUNT,  cudaMemcpyDeviceToHost);
//    if (error_1){
//        std::cout<<"ERROR_11";
//    }
    //фУРЬЕ =============================================
/*
    double after_arr[10000];
    struct d fourier[10000];

    fftw_plan plan;
    plan = fftw_plan_dft_1d(POW2, (fftw_complex *) df->fourfold_phase_data, (fftw_complex *) fourier, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    FILE *file2;
    file2 = fopen("../double_freq.txt", "w");

    for (int i = 0; i < POW2; ++i) {
        after_arr[i] = sqrt(fourier[i].r * fourier[i].r + fourier[i].i * fourier[i].i);

        fprintf(file2, "%lf\n", after_arr[i]);
    }
    fclose(file2);

//    fftw_destroy_plan(plan);
//    ====================================================
*/
}

__global__ void summing_gpu(struct DF * device_df){
    for (int i = 0; i < COUNT; ++i) {
        double tmp_C = 0;
        double tmp_S = 0;
        for (int j = 0; j < 10000; ++j) {
            tmp_C += device_df[i].fourfold_phase_data[j].C;
            tmp_S += device_df[i].fourfold_phase_data[j].S;
        }
        device_df[i].Cs = tmp_C / 10000.0;
        device_df[i].Ss = tmp_S / 10000.0;
        device_df[i].ampl = sqrt(device_df[i].Cs * device_df[i].Cs + device_df[i].Ss * device_df[i].Ss);
        device_df[i].phi = atan2(device_df[i].Ss, device_df[i].Cs);
        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf    %lf\n", i, device_df[i].phi, device_df[i].ampl, device_df[i].Cs, device_df[i].Ss , atan(device_df[i].Ss/device_df[i].Cs));
    }
}

void summing(struct DF * df) {
    for (int i = 0; i < COUNT; ++i) {
        double tmp_C = 0;
        double tmp_S = 0;
        for (int j = 0; j < 10000; ++j) {
            tmp_C += df[i].fourfold_phase_data[j].C;
            tmp_S += df[i].fourfold_phase_data[j].S;
        }
        df[i].Cs = tmp_C / 10000.0;
        df[i].Ss = tmp_S / 10000.0;
        df[i].ampl = sqrt(df[i].Cs * df[i].Cs + df[i].Ss * df[i].Ss);
        df[i].phi = atan2(df[i].Ss, df[i].Cs);
//        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf    %lf\n", i, df[i].phi, df[i].ampl, df[i].Cs, df[i].Ss , atan(df[i].Ss/df[i].Cs));
    }
}

__global__ void multiply(float *in, float *out){
    int tid = blockIdx.x;
    int cur_thread_id;
    for (int i = 0; i < ARR_SIZE; i++){
        cur_thread_id = i % THREADS_GPU_COUNT;
        if (tid < THREADS_GPU_COUNT && cur_thread_id == tid){
            out[i] = in[i] * 2;
        }
    }
}

void multiply_cpu(float *in, float *out){
    for (int i = 0; i < ARR_SIZE; i++){
        out[i] = in[i] * 2;
    }
}

int main()
{
    auto * file_dec = static_cast<Decoder *>(malloc(sizeof(struct Decoder)));
    read_file("../ru0883_sv_no0026.m5b", file_dec);
    for (int i = 0; i < COUNT; ++i) {
        record_float_to_file(file_dec->df[i].decoded_data);
    }

// ЗАГРУЗКА ДАННЫХ В ПАЯМТЬ GPU  ===========================

/*
// ===============ТЕСТ ПРООИЗВОДИТЕЛЬНОСТИ ГПУ =====================
    float test_arr_in[ARR_SIZE];
    float test_arr_out[ARR_SIZE];
    float *dev_test_arr_in;
    float *dev_test_arr_out;
    int error_0 = 0;

    for (int i = 0; i <ARR_SIZE; i++){
        test_arr_in[i] = (float)i;
    }
    error_0 = cudaMalloc((void**) &dev_test_arr_in, ARR_SIZE * sizeof(float));
    if (error_0){
        printf("%d ERROR MEM ALLOCATION\n", error_0);
    }

    error_0 = cudaMalloc((void**) &dev_test_arr_out, ARR_SIZE * sizeof(float));
    if (error_0){
        printf("%d ERROR MEM ALLOCATION\n", error_0);
    }

    error_0 = cudaMemcpy(dev_test_arr_in, test_arr_in, ARR_SIZE * sizeof(float), cudaMemcpyHostToDevice);
    if (error_0){
        printf("%d ERROR MEM COPY TO GPU\n", error_0);
    }
    auto start1 = std::chrono::steady_clock::now();
    multiply<<<THREADS_GPU_COUNT,1>>>(dev_test_arr_in, dev_test_arr_out);
    auto end1 = std::chrono::steady_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    auto start2 = std::chrono::steady_clock::now();
    multiply_cpu(test_arr_in, test_arr_out);
    auto end2 = std::chrono::steady_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);

    error_0 = cudaMemcpy(test_arr_out, dev_test_arr_out, ARR_SIZE * sizeof(float), cudaMemcpyDeviceToHost);
    if (error_0){
        printf("%d ERROR MEM COPY TO CPU\n", error_0);
    }

    for (int i = 0; i < ARR_SIZE; i++){
        printf("%f\n", test_arr_out[i]);
    }

    std::cout << "Trivial mult time GPU: " << (double) duration1.count() << "\n";
    std::cout << "Trivial mult time CPU: " << (double) duration2.count() << "\n";

// =============== ТЕСТ ПРОИЗВОДИТЕЛЬНОСТИ ГПУ ===============
 */


    int error_1 = 0;
    float * data;
    error_1 = cudaMalloc((void**) &data, 10000 * sizeof(float));
    if (error_1){
        printf("%d ERROR MEM ALLOCATION\n", error_1);
    }
//    int error_2 = 0;
//    error_2 = cudaMemcpy(data, file_dec->df[0].decoded_data, 10000 * sizeof(float), cudaMemcpyHostToDevice);
//    if (error_2){
//        printf("%d ERROR MEM COPY TO GPU\n", error_2);
//    }

// ==============================


    phase_rotation(file_dec->df);
//    FilterLowerFreq(file_dec->df);

    struct DF * device_df_lf;
    cudaMalloc((void **)&device_df_lf, sizeof(struct DF) * COUNT);
    cudaMemcpy(device_df_lf, file_dec->df, sizeof(struct DF) * COUNT, cudaMemcpyHostToDevice);
    FilterLowerFreqGPU<<<1,1>>>(device_df_lf);


//    frequency_fourfold(file_dec->df);
    frequency_fourfold_GPU<<<1,1>>>(device_df_lf);
    cudaMemcpy(file_dec->df, device_df_lf, sizeof(struct DF) * COUNT, cudaMemcpyDeviceToHost);
    cudaFree(device_df_lf);
    struct DF * device_df;
    cudaMalloc((void **)&device_df, sizeof(struct DF) * COUNT);
    cudaMemcpy(device_df, file_dec->df, sizeof(struct DF) * COUNT, cudaMemcpyHostToDevice);
    summing_gpu<<<1,1>>>(device_df);
//    summing(file_dec->df);
    cudaMemcpy(file_dec->df, device_df, sizeof(struct DF) * COUNT, cudaMemcpyDeviceToHost);

    free(file_dec);
    cudaFree(device_df);
//    getInfo();
    return EXIT_SUCCESS;
}