#include <stdio.h>
#include "iostream"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "chrono"

//#include "code.cu"

/*#if defined(_WIN32) || defined(_WIN64)
#	include "fftw3.h"
#else

#	include <fftw3.h>

#endif*/

//#define COUNT 4

#define THREADS_GPU_COUNT 5
#define ARR_SIZE 1000000

//#define POW2 8192
#define POW2 10000
#define M_PI 3.14159265358979323846
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
    char data[POW2];
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

/*
double sin_arr[]=
    { 
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
}

*/

/*

*/
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

__global__ void phase_rotation_GPU(struct DF * device_df, float device_cos_arr[], float device_sin_arr[]){

  /*  double alpha = f * (double)t;
        alpha-=floor(alpha);
        alpha*=2*M_PI;
        return alpha;
   */
   
   /* struct Frequency_pair res;
    float f = 1851.0 / 8192.0;
    float alpha=calc_alpha(t,f);
    printf("alpha = %f\n", alpha);
    res.C = I * cos(alpha);
    res.S = I * sin(alpha);
    return res;
*/
    float f = 1851.0 / 8192.0;
    double alpha, alpha_tmp;
    int alpha_int;
    for (int j = 0; j < POW2; ++j) {

        //device_df->phase_data[j].S = calc_freq_pair(j,device_df->decoded_data[j]).S;
        //device_df->phase_data[j].C = calc_freq_pair(j,device_df->decoded_data[j]).C;
        alpha = f * (double)j;
        alpha-=floor(alpha);
        //printf("alpha = %f\n", alpha);
        alpha *=2*M_PI;
        //printf("\n\n%llf\n\n",M_PI);
        alpha_tmp = alpha * 180.0 / M_PI;
        //printf("alpha = %f\n", alpha);
        alpha_int = floor(alpha_tmp);
   //     printf("alpha = %d\n", alpha_int);
    float  cos_arr[360] = { 1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 0.642788, 0.629320, 0.615662, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, 0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069757, 0.052336, 0.034899, 0.017452, -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104529, -0.121869, -0.139173, -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453991, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798635, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948, -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519, -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152, -0.809017, -0.798635, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615662, -0.601815, -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472, -0.453991, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568, -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648, -0.156435, -0.139173, -0.121869, -0.104528, -0.087156, -0.069757, -0.052336, -0.034899, -0.017452, 0.000000, 0.017452, 0.034899, 0.052336, 0.069757, 0.087156, 0.104528, 0.121869, 0.139173, 0.156435, 0.173648, 0.190809, 0.207911, 0.224951, 0.241922, 0.258819, 0.275637, 0.292371, 0.309017, 0.325568, 0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810, 0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320, 0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913546, 0.920505, 0.927184, 0.933580, 0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848};
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

__global__ void FilterLowerFreqGPU(struct DF * device_df){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

//    for(int k = 0; k < COUNT; ++k){
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
                local_avg_re += device_df->phase_data[index].C;
                local_avg_im = device_df->phase_data[index].S;
            }

            device_df->after_filter_data[i].C = local_avg_re / 8.0;
            device_df->after_filter_data[i].S = local_avg_im / 8.0;
        }
    //}


}

void FilterLowerFreq(struct DF * df){
    double local_avg_re=0;
    double local_avg_im=0;
    long int index = 0;

//    for(int k = 0; k < COUNT; ++k){
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
   // }


    //фУРЬЕ =============================================
/*
    double after_arr[POW2];
    struct d fourier[POW2];

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
*/
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

__global__ void frequency_x4_gpu(struct Frequency_pair incoming, struct Frequency_pair * out){
   
    struct Frequency_pair frequency_2;
    frequency_2.C = (incoming.C * incoming.C) - (incoming.S * incoming.S);
    frequency_2.S = 2 * incoming.C * incoming.S;

    out->C = (frequency_2.C * frequency_2.C) - (frequency_2.S * frequency_2.S);
    out->S = 2 * frequency_2.C * frequency_2.S;
}

void frequency_multiply(struct DF * device_df){
//        for (int i = 0; i < COUNT; ++i) {

            for (int j = 0; j < POW2; ++j) {

                frequency_x4_gpu<<<1,1>>>(device_df->after_filter_data[j],&device_df->fourfold_phase_data[j]);
//                df[i].fourfold_phase_data[j].S = frequency_x4_gpu(df[i].after_filter_data[j]).S;
//                df[i].fourfold_phase_data[j].C = frequency_x4_gpu(df[i].after_filter_data[j]).C;
//                df[i].fourfold_phase_data[j].S = frequency_x4_gpu(df[i].after_filter_data[j]).S;
            }
//    }
}

__global__ void frequency_fourfold_GPU(struct DF * device_df){
   // for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < POW2; ++j) {
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
//    struct DF * device_df;
//
//    int error_1;
//    error_1 = cudaMalloc((void**) &device_df, sizeof(struct DF) * COUNT);
//    if (error_1){
//        std::cout<<"ERROR_9";
//    }

   // for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < POW2; ++j) {
            df->fourfold_phase_data[j].C = frequency_x4(df->after_filter_data[j]).C;
            df->fourfold_phase_data[j].S = frequency_x4(df->after_filter_data[j]).S;
        }
    //}
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
   // for (int i = 0; i < COUNT; ++i) {
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
        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf    %lf\n", 1, device_df->phi, device_df->ampl, device_df->Cs, device_df->Ss , atan(device_df->Ss/device_df->Cs));
    //}
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
typedef struct {double re; double im;} complex;

void CFT(int sign, double t, double *xre, double *xim, int npow)
{
    complex cxcs,hold,xa;
    int nmax,i,nn,mm,nw,layer,ii,loc,ll,nw1,j,ij;
    int msk[32];
    double  zz,delta,w;

    for(i = 1, nmax = 1; i <= npow; i++)
        nmax*=2;

    zz = 2.*M_PI*(double)sign/(double)nmax;
    if(sign > 0)
        delta = 1./(t*(double)nmax);
    else
        delta = t;

    msk[1] = nmax/2;
    for(i = 2; i <= npow; i++)
        msk[i]=msk[i-1]/2;

    nn = nmax;
    mm = 2;
    for(layer = 1; layer <= npow; layer++)
    {
        nn = nn/2;
        nw = 0;
        for(i = 1; i<= mm; i+=2)
        {
            ii = nn*i;
            w = zz*(double)nw;
            cxcs.re = cos(w);
            cxcs.im = sin(w);
            for(j = 1; j <= nn; j++)
            {
                ii++;
                ij = ii - nn;
                xa.re = cxcs.re*xre[ii-1] - cxcs.im*xim[ii-1];
                xa.im = cxcs.re*xim[ii-1] + cxcs.im*xre[ii-1];
                xre[ii-1] = xre[ij-1] - xa.re;
                xim[ii-1] = xim[ij-1] - xa.im;
                xre[ij-1] = xre[ij-1] + xa.re;
                xim[ij-1] = xim[ij-1] + xa.im;
            }
            for(loc = 2; loc <= npow; loc++)
            {
                ll = nw - msk[loc];
                if(ll < 0)
                {
                    nw += msk[loc];
                    goto L_40;
                }
                if(ll == 0)
                {
                    nw = msk[loc+1];
                    goto L_40;
                }
                nw = ll;
            }
            L_40:;
        }
        mm *= 2;
    }
    nw = 0;
    for(i = 1; i <= nmax; i++)
    {
        nw1 = nw + 1;
        hold.re = xre[nw1-1];
        hold.im = xim[nw1-1];
        if(nw1-i < 0) goto a;
        if(nw1-i > 0)
        {
            xre[nw1-1] = xre[i-1]*delta;
            xim[nw1-1] = xim[i-1]*delta;
        }
        xre[i-1] = hold.re*delta;
        xim[i-1] = hold.im*delta;
        a:		for(loc = 1; loc <= npow; loc++)
    {
        ll = nw - msk[loc];
        if(ll < 0)
        {
            nw += msk[loc];
            goto L_80;
        }
        if(ll == 0)
        {
            nw = msk[loc+1];
            goto L_80;
        }
        nw = ll;
    }
        L_80:;
    }
}

int main()
{
    auto * file_dec = static_cast<Decoder *>(malloc(sizeof(struct Decoder)));
  

    read_file("./ru0883_bd_no0026.m5b", file_dec);
    //for (int i = 0; i < COUNT; ++i) {
        record_float_to_file(file_dec->df->decoded_data);
    //}
 
    double xre[POW2];
    double xim[POW2];
    //struct d fourier[POW2];

    for (int j = 0; j < POW2; ++j) {

        xre[j] = file_dec->df->decoded_data[j];
        xim[j] = 0;
    }
 
    CFT(-1, 1, xre,xim, 13);
    

    /*fftw_plan plan;

    plan = fftw_plan_dft_1d(8192, (fftw_complex *) fourier, (fftw_complex *) fourier, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(plan);
    fftw_destroy_plan(plan);
*/


    FILE *fp;
    fp = fopen("./results.txt", "w");



    double arr[POW2];
    for (int i = 0; i < POW2; ++i) {
        arr[i] = 0;
    }
    for (int i = 0; i < POW2; ++i) {
        arr[i] = sqrt(xre[i] * xre[i] + xim[i] * xim[i]);
       //arr[i] = sqrt(fourier[i].r * fourier[i].r + fourier[i].i * fourier[i].i);
        fprintf(fp, "%lf\n", arr[i]);
//        printf( "%lf\n", arr[i]);
    }
    fclose(fp);

   

    long int max_num = 0;
    float max_f = arr[0];

    for (int i = 0; i < POW2; ++i) {
            if (arr[i]>max_f){
                max_f = arr[i];
                max_num = i;
            }
        }
    printf("\nfloat = %ld | num = %lf\n", max_num, max_f);
    
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


   // int error_1 = 0;
    //float * data;
    //error_1 = cudaMalloc((void**) &data, 10000 * sizeof(float));
   // if (error_1){
     //   printf("%d ERROR MEM ALLOCATION\n", error_1);
    //}
//    int error_2 = 0;
//    error_2 = cudaMemcpy(data, file_dec->df[0].decoded_data, 10000 * sizeof(float), cudaMemcpyHostToDevice);
//    if (error_2){
//        printf("%d ERROR MEM COPY TO GPU\n", error_2);
//    }

// ==============================

   
//    phase_rotation(file_dec->df);
//    FilterLowerFreq(file_dec->df);

    float  cos_arr[360] = { 1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 0.642788, 0.629320, 0.615662, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, 0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069757, 0.052336, 0.034899, 0.017452, -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104529, -0.121869, -0.139173, -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453991, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798635, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948, -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519, -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152, -0.809017, -0.798635, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615662, -0.601815, -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472, -0.453991, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568, -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648, -0.156435, -0.139173, -0.121869, -0.104528, -0.087156, -0.069757, -0.052336, -0.034899, -0.017452, 0.000000, 0.017452, 0.034899, 0.052336, 0.069757, 0.087156, 0.104528, 0.121869, 0.139173, 0.156435, 0.173648, 0.190809, 0.207911, 0.224951, 0.241922, 0.258819, 0.275637, 0.292371, 0.309017, 0.325568, 0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810, 0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320, 0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913546, 0.920505, 0.927184, 0.933580, 0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848};
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

    struct DF * device_df_lf;
    float  device_cos_arr[360] = { 1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, 0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, 0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, 0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, 0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, 0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, 0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, 0.656059, 0.642788, 0.629320, 0.615662, 0.601815, 0.587785, 0.573576, 0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, 0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, 0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, 0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, 0.104528, 0.087156, 0.069757, 0.052336, 0.034899, 0.017452, -0.000000, -0.017452, -0.034899, -0.052336, -0.069756, -0.087156, -0.104529, -0.121869, -0.139173, -0.156434, -0.173648, -0.190809, -0.207912, -0.224951, -0.241922, -0.258819, -0.275637, -0.292372, -0.309017, -0.325568, -0.342020, -0.358368, -0.374607, -0.390731, -0.406737, -0.422618, -0.438371, -0.453991, -0.469472, -0.484810, -0.500000, -0.515038, -0.529919, -0.544639, -0.559193, -0.573576, -0.587785, -0.601815, -0.615661, -0.629320, -0.642788, -0.656059, -0.669131, -0.681998, -0.694658, -0.707107, -0.719340, -0.731354, -0.743145, -0.754710, -0.766044, -0.777146, -0.788011, -0.798635, -0.809017, -0.819152, -0.829038, -0.838671, -0.848048, -0.857167, -0.866025, -0.874620, -0.882948, -0.891006, -0.898794, -0.906308, -0.913545, -0.920505, -0.927184, -0.933580, -0.939693, -0.945519, -0.951056, -0.956305, -0.961262, -0.965926, -0.970296, -0.974370, -0.978148, -0.981627, -0.984808, -0.987688, -0.990268, -0.992546, -0.994522, -0.996195, -0.997564, -0.998630, -0.999391, -0.999848, -1.000000, -0.999848, -0.999391, -0.998630, -0.997564, -0.996195, -0.994522, -0.992546, -0.990268, -0.987688, -0.984808, -0.981627, -0.978148, -0.974370, -0.970296, -0.965926, -0.961262, -0.956305, -0.951057, -0.945519, -0.939693, -0.933580, -0.927184, -0.920505, -0.913545, -0.906308, -0.898794, -0.891007, -0.882948, -0.874620, -0.866025, -0.857167, -0.848048, -0.838671, -0.829038, -0.819152, -0.809017, -0.798635, -0.788011, -0.777146, -0.766044, -0.754710, -0.743145, -0.731354, -0.719340, -0.707107, -0.694658, -0.681998, -0.669131, -0.656059, -0.642788, -0.629320, -0.615662, -0.601815, -0.587785, -0.573576, -0.559193, -0.544639, -0.529919, -0.515038, -0.500000, -0.484810, -0.469472, -0.453991, -0.438371, -0.422618, -0.406737, -0.390731, -0.374607, -0.358368, -0.342020, -0.325568, -0.309017, -0.292372, -0.275637, -0.258819, -0.241922, -0.224951, -0.207912, -0.190809, -0.173648, -0.156435, -0.139173, -0.121869, -0.104528, -0.087156, -0.069757, -0.052336, -0.034899, -0.017452, 0.000000, 0.017452, 0.034899, 0.052336, 0.069757, 0.087156, 0.104528, 0.121869, 0.139173, 0.156435, 0.173648, 0.190809, 0.207911, 0.224951, 0.241922, 0.258819, 0.275637, 0.292371, 0.309017, 0.325568, 0.342020, 0.358368, 0.374607, 0.390731, 0.406737, 0.422618, 0.438371, 0.453991, 0.469472, 0.484810, 0.500000, 0.515038, 0.529919, 0.544639, 0.559193, 0.573576, 0.587785, 0.601815, 0.615662, 0.629320, 0.642788, 0.656059, 0.669131, 0.681998, 0.694658, 0.707107, 0.719340, 0.731354, 0.743145, 0.754710, 0.766044, 0.777146, 0.788011, 0.798636, 0.809017, 0.819152, 0.829038, 0.838671, 0.848048, 0.857167, 0.866025, 0.874620, 0.882948, 0.891007, 0.898794, 0.906308, 0.913546, 0.920505, 0.927184, 0.933580, 0.939693, 0.945519, 0.951057, 0.956305, 0.961262, 0.965926, 0.970296, 0.974370, 0.978148, 0.981627, 0.984808, 0.987688, 0.990268, 0.992546, 0.994522, 0.996195, 0.997564, 0.998630, 0.999391, 0.999848};
    float  device_sin_arr[360] = { 
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



    cudaMalloc((void **)&cos_arr, 360 * sizeof(float));
    cudaMalloc((void **)&sin_arr, 360 * sizeof(float));

    cudaMemcpy(device_cos_arr, cos_arr, sizeof(cos_arr), cudaMemcpyHostToDevice);
    cudaMemcpy(device_sin_arr, sin_arr, sizeof(sin_arr), cudaMemcpyHostToDevice);



    cudaMalloc((void **)&device_df_lf, sizeof(struct DF));
    cudaMemcpy(device_df_lf, file_dec->df, sizeof(struct DF), cudaMemcpyHostToDevice);

    

    phase_rotation_GPU<<<1,1>>>(device_df_lf, device_cos_arr, device_sin_arr);

    FilterLowerFreqGPU<<<1,1>>>(device_df_lf);
//    frequency_fourfold(file_dec->df);
    frequency_fourfold_GPU<<<1,1>>>(device_df_lf);
    cudaMemcpy(file_dec->df, device_df_lf, sizeof(struct DF), cudaMemcpyDeviceToHost);
    cudaFree(device_df_lf);
    struct DF * device_df;
    cudaMalloc((void **)&device_df, sizeof(struct DF));
    cudaMemcpy(device_df, file_dec->df, sizeof(struct DF), cudaMemcpyHostToDevice);
    summing_gpu<<<1,1>>>(device_df);
//    summing(file_dec->df);
    cudaMemcpy(file_dec->df, device_df, sizeof(struct DF), cudaMemcpyDeviceToHost);

    free(file_dec);
    cudaFree(device_df);
    cudaFree(device_cos_arr);
    cudaFree(device_sin_arr);
//    getInfo();
    return EXIT_SUCCESS;
}