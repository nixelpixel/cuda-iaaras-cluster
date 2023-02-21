#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//#include "code.cu"

#if defined(_WIN32) || defined(_WIN64)
#	include "fftw3.h"
#else

#	include <fftw3.h>

#endif

#define COUNT 4

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

        for (int j = 0; j < 16; ++j) {
//            printf("[%d] = %x\n", j, decoder->df[i].header.hat[j]);
            printf("[%d] data = %x\n", j, decoder->df[i].data[j]);
        }
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


int main()
{
    auto * file_dec = static_cast<Decoder *>(malloc(sizeof(struct Decoder)));
    read_file("../ru0883_bd_no0026.m5b", file_dec);
//    for (int i = 0; i < COUNT; ++i) {
//        record_float_to_file(file_dec->df[i].decoded_data);
////
//    }


    return 0;
}