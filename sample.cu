#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define THREADS_GPU_COUNT 5
#define ARR_SIZE 1000000

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
    //dataunit_t data[POW2];
    float decoded_data[POW2];
    struct Frequency_pair phase_data[POW2];
    struct Frequency_pair after_filter_data[POW2];
    struct Frequency_pair fourfold_phase_data[POW2];
    double Cs, Ss;
    double ampl;
    double phi;
    double f_next;


};

struct DF  p = {};


__global__ void a(struct DF & p){

    p->decoded_data[1] = 4.44;
    printf("\n%f\n",p.decoded_data[1]);
}

__global__ void b(struct DF * p){
    printf("\n%f\n",p.decoded_data[1]);
}

int main(){

   // struct DF * p;

    cudaMalloc((void **) &p, sizeof(struct DF) );

    a<<<1,1>>>(&p);
    b<<<1,1>>>(&p);
    return 0;

}