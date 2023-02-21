#include "CFT.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
#	include "fftw3.h"
#else

#	include <fftw3.h>

#endif

#define WRITE_CMPLX 0 /* what will write in file */

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

void record_float_to_file(float array[10000] ) {

    FILE *fp;
    fp = fopen("../floats.txt", "w");

    for (unsigned i = 0; i < 10000; i++) {
        fprintf(fp, "%f\n", array[i]);

    }

    fclose(fp);

}
//***********************************************************************************************************
// CFT - Complex Fourier Transform and other auxiliary functions
//**********************************************************************************************************


typedef struct {double re;double im;} complex;

//----------------------------------------------------------------------------------------------------------------------

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
//----------------------------------------------------------------------------------------------------------------------

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

void frequency_fourfold(struct DF * df){
    for (int i = 0; i < COUNT; ++i) {
        for (int j = 0; j < 10000; ++j) {
            df[i].fourfold_phase_data[j].C = frequency_x4(df[i].after_filter_data[j]).C;
            df[i].fourfold_phase_data[j].S = frequency_x4(df[i].after_filter_data[j]).S;
        }
    }

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

void print_parametrs(struct DF * df){
    for (int i = 0; i < COUNT; ++i) {
        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf \n");

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
        printf("[%d]: phi = %lf | ampl = %lf | Cs = %lf | Ss = %lf    %lf\n", i, df[i].phi, df[i].ampl, df[i].Cs, df[i].Ss , atan(df[i].Ss/df[i].Cs));
    }
}



int code() {

    struct Decoder * file_dec = static_cast<Decoder *>(malloc(sizeof(struct Decoder)));
    read_file("ru0883_bd_no0026.m5b", file_dec);
    for (int i = 0; i < COUNT; ++i) {
        record_float_to_file(file_dec->df[i].decoded_data);
//
    }


    struct d fft[POW2];
//
//    double xre[POW2];
//    double xim[POW2];


//    for (int i = 0; i < COUNT; ++i) {
    for (int j = 0; j < POW2; ++j) {

        fft[j].r = file_dec->df[0].decoded_data[j];
        fft[j].i = 0;
//            xre[i] = file_dec->df[0].decoded_data[i];
//            xim[i] = 0;
    }
//    }

    // ДО ФАЗОВРАЩЕНИЯ
    /*   fftw_plan plan;
       plan = fftw_plan_dft_1d(POW2, (fftw_complex *) fft, (fftw_complex *) fft, FFTW_FORWARD, FFTW_ESTIMATE);
       fftw_execute(plan);
   */
//    CFT(1,0.001,xre,xim,10);

    FILE *fp;
    fp = fopen("../results.txt", "w");

//    for (unsigned i = 0; i < 10; i++) {
//        printf( "%lf %lf\n", fft[i].r, fft[i].i);
//    }

    double arr[POW2];
    for (int i = 0; i < POW2; ++i) {
        arr[i] = 0;
    }

    long int max_num = 0;
    float max_f = arr[0];



    for (int i = 0; i < POW2; ++i) {
        arr[i] = sqrt(fft[i].r * fft[i].r + fft[i].i * fft[i].i);

        fprintf(fp, "%lf\n", arr[i]);
//        printf( "%lf\n", arr[i]);
    }

    fclose(fp);
//    fftw_destroy_plan(plan);

    phase_rotation(file_dec->df);

    for (unsigned i = 0; i < 10; i++) {
        printf( "%lf %lf\n", file_dec->df[0].phase_data[i].C, file_dec->df[0].phase_data[i].S);
    }
    // ПОСЛЕ ФАЗАВРАЩЕНИЯ
/*
//    plan = fftw_plan_dft_1d(POW2, (fftw_complex *) file_dec->df->phase_data, (fftw_complex *) file_dec->df->phase_data, FFTW_FORWARD, FFTW_ESTIMATE);
//    fftw_execute(plan);
//
//    FILE *file;
//    file = fopen("../after_phase_shift.txt", "w");
//
//    double phase_arr[POW2];
//    for (int i = 0; i < POW2; ++i){
//        phase_arr[i] = sqrt(file_dec->df[0].phase_data[i].C * file_dec->df[0].phase_data[i].C + file_dec->df[0].phase_data[i].S * file_dec->df[0].phase_data[i].S);
//        fprintf(file, "%lf\n", phase_arr[i]);
//    }
//    fclose(file);
*/

    FilterLowerFreq(file_dec->df);
    frequency_fourfold(file_dec->df);
    summing(file_dec->df);

    for (int i = 0; i < POW2; ++i) {
        if (arr[i]>max_f){
            max_f = arr[i];
            max_num = i;
        }
    }
    printf("\nfloat = %ld | num = %lf\n", max_num, max_f);


//    fftw_destroy_plan(plan);

    free(file_dec);
}
