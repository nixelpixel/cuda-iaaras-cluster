#include <stdio.h>
//#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "libraries/cuda_code.h"
#include "libraries/zclog_lib.h"
#include "libraries/general_types.h"
#include <time.h>


//#define COUNT 4

#define THREADS_GPU 5
#define BLOCKS_GPU 10000

#define ARR_SIZE 1000000

//#define POW2 8192
#define POW2 10000
#define M_PI 3.14159265358979323846

#ifdef __cplusplus
//extern "C" {
#endif

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
    dataunit_t data[POW2];
    float decoded_data[POW2];
    struct Frequency_pair phase_data[POW2];
    struct Frequency_pair after_filter_data[POW2];
    struct Frequency_pair fourfold_phase_data[POW2];
    double Cs, Ss;
    double ampl;
    double phi;
    double f_next;


};

struct Decoder {

    FILE *file;
    struct DF * df;

};

struct DF * device_memory_pointer;

float * sinus_array;
float * cosinus_array;

void getInfo(){
    int count;
    cudaGetDeviceCount(&count);
    cudaDeviceProp prop = {};
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

const float const_DecodeAmlitudeArray[]={
    -0.3333, -1.0, 0.3333, 1.0
};


void decoding(char in, float * res){
    unsigned int amplindex = ( in >> 2 ) & 0x03;

    *res = const_DecodeAmlitudeArray[ amplindex ];
}
/*
void decoding(char in, float * res){

    char arr;

    arr = (in >> 2) & 0x3;

    if (arr == 0x0){
        *res = -0.3333;
    }
    if (arr == 0x1){
        *res = -1.0;
    }
    if (arr == 0x2){
        *res = 0.3333;
    }
    if (arr == 0x3){
        *res = 1.0;
    }

}*/

void read_file(char *file_name, struct Decoder *decoder) {

    decoder->file = fopen(file_name, "rb");
    if (!decoder->file) {
        perror("File error");
        return;
    }

    FILE *f = decoder->file;
    decoder->df = (struct DF *) malloc(sizeof(struct DF));
    //for (int i = 0; i < COUNT; ++i) {
        fread(decoder->df->header.hat, 16, 1, f);
        fread(decoder->df->data, POW2, 1, f);

        for (int j = 0; j < 16; ++j) {
            printf("[%d] = %x\n", j, (unsigned char) decoder->df->header.hat[j]);
//            printf("[%d] data = %x\n", j, decoder->df[i].data[j]);
        }
        printf("\n");
    //}

    //DECODING
   // for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < POW2; ++j) {
            decoding(decoder->df->data[j],&decoder->df->decoded_data[j]);
//            printf("data[%d] = %f\n",i, decoder->df[j].decoded_data[i]);
        }
        fclose(decoder->file);
    //}
}




void record_float_to_file(float array[POW2] ) {

    FILE *fp;
    fp = fopen("../floats.txt", "w");

    for (unsigned i = 0; i < POW2; i++) {
        fprintf(fp, "%f\n", array[i]);

    }

    fclose(fp);

}

float calc_alpha(long int t, float f){

    double alpha = f * (double)t;
    alpha-=floor(alpha);
    alpha*=2* M_PI;
    return alpha;
}

struct Frequency_pair calc_freq_pair(long int t, float I){
    struct Frequency_pair res;
    float f = 1851.0 / 8192.0;
    float alpha=calc_alpha(t,f);
    printf("alpha = %f\n", alpha);
    res.C = I * cos(alpha);
    res.S = I * sin(alpha);
    return res;
}

/*__global__ struct Frequency_pair calc_freq_pair_GPU(long int t, float I){
    struct Frequency_pair res;
    float f = 1851.0 / 8192.0;
    float alpha=calc_alpha(t,f);
    printf("alpha = %f\n", alpha);
    res.C = I * cos(alpha);
    res.S = I * sin(alpha);
    return res;
}*/

void phase_rotation(struct DF * df) {
    FILE *fp;
    fp = fopen("../phase.txt", "w");
    //for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < POW2; ++j) {
            df->phase_data[j].S = calc_freq_pair(j,df->decoded_data[j]).S;
            df->phase_data[j].C = calc_freq_pair(j,df->decoded_data[j]).C;
//            printf("phase[%d][%d]: C = %f| S = %f\n", i,j,df[i].phase_data[j].C,df[i].phase_data[j].S);

            fprintf(fp, "%f %f\n",df->phase_data[j].C,df->phase_data[j].S);

        }

   // }
    fclose(fp);

}

__global__ void phase_rotation_GPU(struct DF * device_df, double freq, datablockid_t datablock_id/*, float * device_cos_arr, float * device_sin_arr*/){

    double alpha, alpha_tmp;
    int alpha_int;
    int j = blockIdx.x;
      
    if (j < BLOCKS_GPU){
    //for (int j = 0; j < POW2; ++j) {

        alpha = freq * (double)j;
        alpha-=floor(alpha);

        alpha *=2*M_PI;
   
        alpha_tmp = alpha * 180.0 / M_PI;

        alpha_int = floor(alpha_tmp);

        
        float  cos_arr[360] = { 
            1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688,
            0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519,
            0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620,
            0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146,
            0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059,
            0.642788, 0.629320, 0.615662, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038,
            0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368,
            0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809,
            0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069757, 0.052336, 0.034899, 0.017452,
            -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104529, -0.121869, -0.139173,
            -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372,
            -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371,
            -0.453991, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576,
            -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658,
            -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798635,
            -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948,
            -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519,
            -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808,
            -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848,
            -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268,
            -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305,
            -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794,
            -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152,
            -0.809017, -0.798635, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340,
            -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615662, -0.601815,
            -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472,
            -0.453991, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568,
            -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648,
            -0.156435, -0.139173, -0.121869, -0.104528, -0.087156, -0.069757, -0.052336, -0.034899, -0.017452,
            0.000000, 0.017452, 0.034899, 0.052336, 0.069757, 0.087156, 0.104528, 0.121869, 0.139173, 0.156435,
            0.173648, 0.190809, 0.207911, 0.224951, 0.241922, 0.258819, 0.275637, 0.292371, 0.309017, 0.325568,
            0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810,
            0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320,
            0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710,
            0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167,
            0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913546, 0.920505, 0.927184, 0.933580,
            0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627,
            0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848
        };

        float  sin_arr[360] = { 
            0.000000, 0.017452, 0.034899, 0.052336, 0.069756, 0.087156, 0.104528, 0.121869, 0.139173, 0.156434, 
            0.173648, 0.190809, 0.207912, 0.224951, 0.241922, 0.258819, 0.275637, 0.292372, 0.309017, 0.325568, 
            0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810, 
            0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320, 
            0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 
            0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 
            0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913545, 0.920505, 0.927184, 0.933580, 
            0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 
            0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848, 
            1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 
            0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951056, 0.945519, 
            0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891006, 0.882948, 0.874620, 
            0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798635, 0.788011, 0.777146, 
            0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 
            0.642788, 0.629321, 0.615661, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 
            0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374606, 0.358368, 
            0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 
            0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069756, 0.052336, 0.034899, 0.017452, 
            -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104528, -0.121869, -0.139173, -0.156434, 
            -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, 
            -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453991, -0.469472, -0.484810, 
            -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, 
            -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, 
            -0.766045, -0.777146, -0.788011, -0.798635, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, 
            -0.866025, -0.874620, -0.882948, -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933581, 
            -0.939693, -0.945519, -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, 
            -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, 
            -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, 
            -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951056, -0.945519, 
            -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891006, -0.882948, -0.874620, 
            -0.866025, -0.857167, -0.848048, -0.838670, -0.829038, -0.819152, -0.809017, -0.798635, -0.788011, -0.777146, 
            -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, 
            -0.642788, -0.629320, -0.615661, -0.601815, -0.587785, -0.573577, -0.559193, -0.544639, -0.529919, -0.515038, 
            -0.500000, -0.484809, -0.469471, -0.453991, -0.438371, -0.422618, -0.406736, -0.390731, -0.374607, -0.358368, 
            -0.342020, -0.325568, -0.309017, -0.292372, -0.275638, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, 
            -0.173648, -0.156434, -0.139173, -0.121869, -0.104529, -0.087156, -0.069756, -0.052336, -0.034900, -0.017453
        };
        
       // printf("cos = %f\n", cos_arr[alpha_int]);
        device_df->phase_data[j].C = device_df->decoded_data[j] * cos_arr[alpha_int];
        device_df->phase_data[j].S = device_df->decoded_data[j] * sin_arr[alpha_int];
        

    }

    
}


__global__ void FilterLowerFreq_GPU(struct DF * device_df, datablockid_t datablock_id){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

    //printf("F[%u]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf\n", datablock_id, device_df->phi, device_df->ampl, device_df->Cs, device_df->Ss);

    //for (int i = 0; i < POW2; ++i) {
    int i = blockIdx.x;
      
    if (i < BLOCKS_GPU){
        local_avg_re = 0;
        local_avg_im = 0;

        for (int j = -4; j <= 3; ++j) {
            index = i+j;
            if (index < 0){
                index = 0;
            }
            else if ( index >= POW2) {
                index = 9999;
            }
            local_avg_re += device_df->phase_data[index].C;
            local_avg_im = device_df->phase_data[index].S;
        }

        device_df->after_filter_data[i].C = local_avg_re / 8.0;
        device_df->after_filter_data[i].S = local_avg_im / 8.0;
    }

}

void FilterLowerFreq(struct DF * df){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

    for (int i = 0; i < POW2; ++i) {
        local_avg_re = 0;
        local_avg_im = 0;

        for (int j = -4; j <= 3; ++j) {
            index = i+j;
            if (index < 0){
                index = 0;
            }
            else if ( index >= POW2) {
                index = 9999;
            }
            local_avg_re += df->phase_data[index].C;
            local_avg_im = df->phase_data[index].S;
        }

        df->after_filter_data[i].C = local_avg_re / 8.0;
        df->after_filter_data[i].S = local_avg_im / 8.0;
    }
}

struct Frequency_pair frequency_x4(struct Frequency_pair incoming) {

    struct Frequency_pair frequency_2;
    frequency_2.C = (incoming.C * incoming.C) - (incoming.S * incoming.S);
    frequency_2.S = 2 * incoming.C * incoming.S;

    struct Frequency_pair frequency_4;
    frequency_4.C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
    frequency_4.S = 2 * frequency_2.C * frequency_2.S;
    return frequency_4;
}

// __global__ void frequency_x4_gpu(struct Frequency_pair incoming, struct Frequency_pair * out){
   
//     struct Frequency_pair frequency_2;
//     frequency_2.C = (incoming.C * incoming.C) - (incoming.S * incoming.S);
//     frequency_2.S = 2 * incoming.C * incoming.S;

//     out->C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
//     out->S = 2 * frequency_2.C * frequency_2.S;
// }

// void frequency_multiply(struct DF * device_df){

//     for (int j = 0; j < POW2; ++j) {
//         frequency_x4_gpu<<<1,1>>>(device_df->after_filter_data[j],&device_df->fourfold_phase_data[j]);
//     }
// }

__global__ void frequency_fourfold_GPU(struct DF * device_df){
    
    int j = blockIdx.x;

    //for (int j = 0; j < BLOCKS_GPU; ++j) {
    if (j < BLOCKS_GPU){
        struct Frequency_pair frequency_2;
        frequency_2.C = (device_df->after_filter_data[j].C * device_df->after_filter_data[j].C) - (device_df->after_filter_data[j].S * device_df->after_filter_data[j].S);
        frequency_2.S = 2 * device_df->after_filter_data[j].C * device_df->after_filter_data[j].S;

        struct Frequency_pair frequency_4;
        frequency_4.C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
        frequency_4.S = 2 * frequency_2.C * frequency_2.S;

        device_df->fourfold_phase_data[j].C = frequency_4.C;
        device_df->fourfold_phase_data[j].S = frequency_4.S;
    }
    //}
}

void frequency_fourfold(struct DF * df){


        for (int j = 0; j < POW2; ++j) {
            df->fourfold_phase_data[j].C = frequency_x4(df->after_filter_data[j]).C;
            df->fourfold_phase_data[j].S = frequency_x4(df->after_filter_data[j]).S;
        }
}

__global__ void summing_gpu(struct DF * device_df, datablockid_t datablock_id){
	//zclog( (char *) "$sizeof( struct DF ): %i bytes.", sizeof( device_df ) );

  
        double tmp_C = 0;
        double tmp_S = 0;
        for (int j = 0; j < POW2; ++j) {
            tmp_C += device_df->fourfold_phase_data[j].C;
            tmp_S += device_df->fourfold_phase_data[j].S;
        }
        device_df->Cs = tmp_C / POW2;
        device_df->Ss = tmp_S / POW2;
        device_df->ampl = sqrt(device_df->Cs * device_df->Cs + device_df->Ss * device_df->Ss);
        device_df->phi = atan2(device_df->Ss, device_df->Cs);
        //printf("[%u]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf\n", datablock_id, device_df->phi, device_df->ampl, device_df->Cs, device_df->Ss);
	printf("[%u]: phi = %lf | ampl = %lf\n", datablock_id, (device_df->phi)*180.0/M_PI, device_df->ampl);
   	

}

void summing(struct DF * df) {
    //for (int i = 0; i < COUNT; ++i) {
        double tmp_C = 0;
        double tmp_S = 0;
        for (int j = 0; j < POW2; ++j) {
            tmp_C += df->fourfold_phase_data[j].C;
            tmp_S += df->fourfold_phase_data[j].S;
        }
        df->Cs = tmp_C / POW2;
        df->Ss = tmp_S / POW2;
        df->ampl = sqrt(df->Cs * df->Cs + df->Ss * df->Ss);
        df->phi = atan2(df->Ss, df->Cs);
//        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf    %lf\n", i, df[i].phi, df[i].ampl, df[i].Cs, df[i].Ss , atan(df[i].Ss/df[i].Cs));
    //}
}

// __global__ void multiply(float *in, float *out){
//     int tid = blockIdx.x;
//     int cur_thread_id;
//     for (int i = 0; i < ARR_SIZE; i++){
//         cur_thread_id = i % THREADS_GPU_COUNT;
//         if (tid < THREADS_GPU_COUNT && cur_thread_id == tid){
//             out[i] = in[i] * 2;
//         }
//     }
// }

// void multiply_cpu(float *in, float *out){
//     for (int i = 0; i < ARR_SIZE; i++){
//         out[i] = in[i] * 2;
//     }
// }
typedef struct {double re; double im;} complex;

struct DF * get_pointer(){
    return device_memory_pointer;
}
/*
__global__ void initializing_triganometry_arrays(float * s_arr, float * c_arr){
    tmp_s_arr = { 
        0.000000, 0.017452, 0.034899, 0.052336, 0.069756, 0.087156, 0.104528, 0.121869, 0.139173, 0.156434, 
        0.173648, 0.190809, 0.207912, 0.224951, 0.241922, 0.258819, 0.275637, 0.292372, 0.309017, 0.325568, 
        0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810, 
        0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320, 
        0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 
        0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 
        0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913545, 0.920505, 0.927184, 0.933580, 
        0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 
        0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848, 
        1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 
        0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951056, 0.945519, 
        0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891006, 0.882948, 0.874620, 
        0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798635, 0.788011, 0.777146, 
        0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 
        0.642788, 0.629321, 0.615661, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 
        0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374606, 0.358368, 
        0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 
        0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069756, 0.052336, 0.034899, 0.017452, 
        -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104528, -0.121869, -0.139173, -0.156434, 
        -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, 
        -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453991, -0.469472, -0.484810, 
        -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, 
        -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, 
        -0.766045, -0.777146, -0.788011, -0.798635, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, 
        -0.866025, -0.874620, -0.882948, -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933581, 
        -0.939693, -0.945519, -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, 
        -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, 
        -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, 
        -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951056, -0.945519, 
        -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891006, -0.882948, -0.874620, 
        -0.866025, -0.857167, -0.848048, -0.838670, -0.829038, -0.819152, -0.809017, -0.798635, -0.788011, -0.777146, 
        -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, 
        -0.642788, -0.629320, -0.615661, -0.601815, -0.587785, -0.573577, -0.559193, -0.544639, -0.529919, -0.515038, 
        -0.500000, -0.484809, -0.469471, -0.453991, -0.438371, -0.422618, -0.406736, -0.390731, -0.374607, -0.358368, 
        -0.342020, -0.325568, -0.309017, -0.292372, -0.275638, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, 
        -0.173648, -0.156434, -0.139173, -0.121869, -0.104529, -0.087156, -0.069756, -0.052336, -0.034900, -0.017453
};
    tmp_c_arr = { 
    1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688,
    0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519,
    0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620,
    0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146,
    0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059,
    0.642788, 0.629320, 0.615662, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038,
    0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368,
    0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809,
    0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069757, 0.052336, 0.034899, 0.017452,
    -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104529, -0.121869, -0.139173,
    -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372,
    -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371,
    -0.453991, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576,
    -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658,
    -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798635,
    -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948,
    -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519,
    -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808,
    -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848,
    -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268,
    -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305,
    -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794,
    -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152,
    -0.809017, -0.798635, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340,
    -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615662, -0.601815,
    -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472,
    -0.453991, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568,
    -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648,
    -0.156435, -0.139173, -0.121869, -0.104528, -0.087156, -0.069757, -0.052336, -0.034899, -0.017452,
    0.000000, 0.017452, 0.034899, 0.052336, 0.069757, 0.087156, 0.104528, 0.121869, 0.139173, 0.156435,
    0.173648, 0.190809, 0.207911, 0.224951, 0.241922, 0.258819, 0.275637, 0.292371, 0.309017, 0.325568,
    0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810,
    0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320,
    0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710,
    0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167,
    0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913546, 0.920505, 0.927184, 0.933580,
    0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627,
    0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848
};
    cudaMemcpy(s_arr, tmp_s_arr, 360 * sizeof( float ), cudaMemcpyHostToDevice);
    cudaMemcpy(c_arr, tmp_c_arr, 360 * sizeof( float ), cudaMemcpyHostToDevice);
}*/

__global__ void initializing_data_gpu(struct DF * device_df /*,float * s_arr, float * c_arr*/){
    device_df->Cs = 0.0;
    device_df->Ss = 0.0;
    device_df->ampl = 0.0;
    device_df->phi = 0.0;
    device_df->f_next = 0.0;
}

clustererrorcode_t gpuglobal_init(){
	return CERR_NONE;
}


clustererrorcode_t gpu_init( corrid_t corrindex, char *logfilename ){
    
    int error_number = 0;
    int count;
    error_number = cudaGetDeviceCount(&count);

    //necessary for "double" calculate
    cudaDeviceProp ness_prop;
    int dev;
    memset( &ness_prop, 0, sizeof( cudaDeviceProp ));
    ness_prop.major = 1;
    ness_prop.minor = 3;
    error_number = cudaChooseDevice( &dev, &ness_prop );
    if(error_number){
        printf("\nCant find requiered device\n");
        return CERR_CUDA;
    }
    dev = 0; // DEVICE 0 or 1
    error_number = cudaSetDevice( dev );
    if(error_number){
        printf("\nCant set device\n");
        return CERR_CUDA;
    }
    else { 
        printf("\nSetet device: %d for correlator: %u\n",dev, corrindex);
    }
    
    error_number = cudaMalloc((void**)&device_memory_pointer, sizeof(struct DF));
    if(error_number){
        printf("\nCant allocate memory\n");
        return CERR_CUDA;
    }
    error_number = cudaMalloc((void**)&sinus_array, 360 * sizeof(float));
    if(error_number){
        printf("\nCant allocate memory\n");
        return CERR_CUDA;
    }

    error_number = cudaMalloc((void**)&cosinus_array, 360 * sizeof(float));
    if(error_number){
        printf("\nCant allocate memory\n");
        return CERR_CUDA;
    }

    //device_memory_pointer->ampl = 0;
    //cudaMemcpy(device_memory_pointer->decoded_data, tempdecoded, 10000, cudaMemcpyHostToDevice);

    initializing_data_gpu<<<1,1>>>(device_memory_pointer);
    //initializing_triganometry_arrays<<<1,1>>>( sinus_array, cosinus_array);
   
    return CERR_NONE;
    
}



double gpu_get_phi( double freq, dataunit_t *data,  datablockid_t datablock_id ) {

    float tempdecoded[ POW2 ];
    for (int j = 0; j < 10000; ++j) {
        //printf( "j: %i. device_df_lf %p, ->data %p, ->decod %p\n", j, device_df_lf, device_df_lf->data, device_df_lf->decoded_data );
        decoding( data[ j ], &tempdecoded[ j ] );
    }
    cudaMemcpy(device_memory_pointer->decoded_data, tempdecoded, 10000 * sizeof( float ), cudaMemcpyHostToDevice);

    double phi;
    int msec = 0, trigger = 10; /* 10ms */
clock_t before = clock();
    phase_rotation_GPU<<<BLOCKS_GPU,1>>>(device_memory_pointer, freq, datablock_id);
    FilterLowerFreq_GPU<<<1,1>>>(device_memory_pointer, datablock_id);
    frequency_fourfold_GPU<<<BLOCKS_GPU,1>>>(device_memory_pointer);
    summing_gpu<<<1,1>>>(device_memory_pointer, datablock_id);
    cudaMemcpy(&phi, &device_memory_pointer->phi, sizeof( device_memory_pointer->phi ) , cudaMemcpyDeviceToHost);
clock_t difference = clock() - before;
  msec = difference * 1000 / CLOCKS_PER_SEC;
printf("Time taken %d seconds %d milliseconds\n",
  msec/1000, msec%1000);    

return phi;

}

struct DF * get_pointer(struct DF * device_df){
    return device_df;
}

void gpuglobal_finalize(){
    printf("\nFinalize\n");
    cudaFree(device_memory_pointer);
    printf("\nMemory has been free...\n");
}

void gpu_finalize( corrid_t corrindex ){}

#ifdef __cplusplus
//}
#endif

/*int main()
{
    struct Decoder *file_dec = (struct Decoder *)malloc(sizeof(struct Decoder));
    read_file("ru0883_bd_no0026.m5b", file_dec);
    get_phi(file_dec->df->data);

    return 0;
}*/
