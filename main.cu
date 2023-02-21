#include "stdio.h"
#include "code.cu"
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

int main()
{
//    int a,b,c;
//    int *dev_c;
//    a=3;
//    b=4;
//    int er = cudaMalloc((void**)&dev_c, sizeof(int));
//    add<<<1,1>>>(a,b,dev_c);
//    cudaMemcpy(&c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);
//    printf("%d + %d is %d\n", a, b, c);
//    printf("error = %d\n", er);
//    cudaFree(dev_c);

//    printf("\ncount =%d \nerror count = %d\n\n", count,error);

//    getInfo();



    return 0;
}