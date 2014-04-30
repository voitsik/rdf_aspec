/***************************************************************************
 *   Copyright (C) 2011-2013 by Petr Voytsik                               *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************/


#define _FILE_OFFSET_BITS 64
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include <complex.h>
#include <fftw3.h>

#define handle_error(msg) \
        do { perror(msg); exit(EXIT_FAILURE); } while (0)


const char program[] = "rdf_aspec";
const char author[]  = "Petr Voytsik";
const char version[] = "2.0";


/* Decoded RDF-file header */
typedef struct rdf_header{
    char sig[5];          /* RDF1 */
    uint16_t header_size; /* Should be 256 bytes */
    time_t date;          /* Unix time of recording start */
    char station[12];
    char source[12];
    char exper[12];
    int data_rate;
    char rec_mode;
    char rdr_mode[4];
} rdf_header_t;

static float lut1bit[256][8];
static float lut2bit[256][4];

static void initluts()
{
	unsigned b, i, l, s, m;
    const float HiMag = 3.3359;
	const float lut2level[2] = {-1.0, 1.0};
    const float lut4level[4] = {-HiMag, 1.0, -1.0, HiMag};

	for(b = 0; b < 256; ++b){
		/* lut1bit */
		for(i = 0; i < 8; ++i){
			l = (b>>i)&1;
			lut1bit[b][i] = lut2level[l];
		}

		/* lut2bit */
		for(i = 0; i < 4; i++){
		    s = i*2;    /* 0, 2, 4, 6 */
		    m = s+1;    /* 1, 3, 5, 7 */
		    l = ((b>>s)&1) + (((b>>m)&1)<<1);
		    lut2bit[b][i] = lut4level[l];
		}
	}
}

/**
 *  decode_1bit
 *  Convert 1bit encoded raw data to 4 float arrays
 *  n - number of samples to decode
 */
static size_t decode_1bit(const uint8_t * const in, float **out, size_t n)
{
    size_t  i;
    float *fp;

    for(i = 0; i < n/2; ++i){
        fp = lut1bit[in[i]];
        out[0][2*i]   = fp[0];
        out[1][2*i]   = fp[1];
        out[2][2*i]   = fp[2];
        out[3][2*i]   = fp[3];
        out[0][2*i+1] = fp[4];
        out[1][2*i+1] = fp[5];
        out[2][2*i+1] = fp[6];
        out[3][2*i+1] = fp[7];
    }

    return n;
}

/**
 *  decode_2bit
 *  Convert 2bit encoded raw data to 4 float arrays
 *  n - number of samples to decode
 */
static size_t decode_2bit(const uint8_t * const in, float **out, size_t n)
{
    size_t  i;
    float *fp;

    for(i = 0; i < n; ++i){
        fp = lut2bit[in[i]];
        out[0][i]   = fp[0];
        out[1][i]   = fp[1];
        out[2][i]   = fp[2];
        out[3][i]   = fp[3];
    }

    return n;
}

#define GET_FIELD(dest, src, n)     \
        do {                        \
            strncpy(dest, src, n);  \
            dest[n] = 0;            \
        } while(0)

/**
 *  parse_rdf_header 
 *  Parse RDF-file header
 *  Check the header signature and size
 *
 *  On success, parse_rdf_header returns the number of bits per sample
 *  On error, -1 is returned
 */
static int parse_rdf_header(const char *h, rdf_header_t *info)
{
    char data_date[19], datarate[3];
    struct tm tm0;
    int bits = 0;

    GET_FIELD(info->sig, h, 4);

    if(strcmp(info->sig, "RDF1")){
        fprintf(stderr, "Wrong RDF signature '%s' while RDF1 expected.\n", info->sig);

        return -1;
    }

    memcpy(&info->header_size, &h[4], 2);
    if(info->header_size != 256){
        fprintf(stderr, "Header size is %u bytes while 256 expected.\n",
                        info->header_size);

        return -1;
    }

    GET_FIELD(data_date, &h[6], 18);

    /* Decode date/time */
    memset(&tm0, 0, sizeof(tm0));
    strptime(data_date, "%Y %j-%H:%M:%S", &tm0);
    info->date = mktime(&tm0);

    GET_FIELD(info->station, &h[24], 11);
    GET_FIELD(info->source, &h[35], 11);
    GET_FIELD(info->exper, &h[46], 11);
    GET_FIELD(datarate, &h[57], 2);
    info->data_rate = atoi(datarate);
    info->rec_mode = h[59];
    GET_FIELD(info->rdr_mode, &h[60], 3);

    fprintf(stderr, " station = %s\n", info->station);
    fprintf(stderr, " source = %s\n", info->source);
    fprintf(stderr, " experiment = %s\n", info->exper);
    fprintf(stderr, " data rate = %d\n", info->data_rate);

    if(info->data_rate == 16 && !strcmp(info->rdr_mode, "DEC")){
        bits = 1;
    }else if(info->data_rate == 32){
        bits = 2;
    }else{
        fprintf(stderr, "Unsupported combination of RDR mode (%s) and data rate %d\n",
                        info->rdr_mode, info->data_rate);
        bits = -1;
    }
    
    fprintf(stderr, " bits = %d\n", bits);

    return bits;
}

/**
 *  print_header prints experiment name and time to stdout
 */
static void print_header(rdf_header_t *h, time_t t_off)
{
    char time_str[64];
    time_t t = h->date + t_off;

    strftime(time_str, sizeof(time_str), "%F %T", gmtime(&t));
    printf("#%s (%s UT)\n", h->exper, time_str);
}

static void usage()
{
    fprintf(stderr, "Usage: %s <RDF-file> <nchan> <nint> [offset]\n", program);
    fprintf(stderr, "  <RDF-file> \t name of input file\n");
    fprintf(stderr, "  <nchan> \t number of channels to make per IF\n");
    fprintf(stderr, "  <nint> \t number of spectra to integrate\n");
    fprintf(stderr, "  [offset] \t number of seconds into file to start decoding\n");
}

int main(int argc, char *argv[])
{
    FILE *f;
    size_t N, M, n;
    off_t offset = 0;
    uint_fast32_t i, j, k;
    uint8_t *raw_data;       /* Raw input data from RDF-file */
    fftwf_complex *c_data[4];  /* Output complex spectrum */
    float *f_data[4];        /* Decoded input data (float) */
    float *sp[4];            /* Power spectrum */
    fftwf_plan p;
    float re, im;
    double df, freq, time_off = 0.;
    int bits;   /* Number of bits per sample */
    rdf_header_t rdf_info;
    char header[256];
    size_t (*decode)(const uint8_t * const , float **, size_t );
    size_t count = 0;
    
    if(argc < 4){
        usage();

        return EXIT_FAILURE;
    }

    f = fopen(argv[1], "r");
    if(!f){
        fprintf(stderr, "Could not open file '%s': ", argv[1]);
        handle_error("");
    }

    /* Number of spectral channels */
    N = atol(argv[2]);
    /* Number of spectra to summ */
    M = atol(argv[3]);

    if(argc == 5){
        time_off = atof(argv[4]);
        offset = (off_t)(400. * 40000. * time_off);
    }

    n = fread(header, sizeof(uint8_t), 256UL, f);
    if(n != 256){
        fprintf(stderr, "Read only %lu of %lu bytes from RDF-file header. Exiting.\n", 
                n, 256UL);
        fclose(f);
        exit(EXIT_FAILURE);
    }

    bits = parse_rdf_header(header, &rdf_info);
    if(bits == 1){
        decode = &decode_1bit;
    }else if(bits == 2){
        decode = &decode_2bit;
    }else{
        if(bits < 0)
            fprintf(stderr, "Error parsing RDF-file header. \n");
        else
            fprintf(stderr, "Wrong number of bits %d\n", bits);

        fclose(f);

        exit(EXIT_FAILURE);
    }

    if(offset)
        if(fseeko(f, offset, SEEK_CUR) < 0)
            handle_error("fseeko");

    raw_data = (uint8_t *)malloc(sizeof(uint8_t)*bits*N);

    for(j = 0; j < 4; j++){
        sp[j] = calloc(N+1, sizeof(float));
        /*c_data[j] = fftwf_alloc_complex(N/2+1);*/ /* Available only since FFTW 3.3-beta1 */
        c_data[j] = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (N+1));
        f_data[j] = (float *)c_data[j];
    }

    if(fftwf_import_system_wisdom() != 1)
        fprintf(stderr, "Warning: could not load system fftw wisdom\n");

    p = fftwf_plan_dft_r2c_1d(2*N, f_data[0], c_data[0], FFTW_ESTIMATE);

    initluts();
    print_header(&rdf_info, time_off);

    for(k = 0; k < M; k++){
        n = fread(raw_data, sizeof(uint8_t), bits*N, f);
        if(n != bits*N){
            perror("fread");
            break;
        }
        
        count += decode(raw_data, f_data, 2*N);

        /* Four channels */
        for(j = 0; j < 4; j++){
            fftwf_execute_dft_r2c(p, f_data[j], c_data[j]);
            for(i = 0; i <= N; i++){
                re = creal(c_data[j][i]);
                im = cimag(c_data[j][i]);
                sp[j][i] += (re*re + im*im)/(float)(2*N);
            }
        }

        if(k % 100 == 0)
            fprintf(stderr, "%.1f%%\r", 100.f * (float)k / (float)M);
    }

    fprintf(stderr, "%lu samples decoded\n", count);
    
    free(raw_data);
    fclose(f);

    /* Clean fftw structures */
    fftwf_destroy_plan(p);
    for(j = 0; j < 4; j++){
        fftwf_free(c_data[j]);
    }
    fftwf_cleanup();

    freq = 0.;
    df = 1. / (1e6 * 31.25e-9 * (double)(2*N));

    for(i = 0; i <= N; i++){
        printf("%f ", freq);
        for(j = 0; j < 4; j++)
            printf("%f ", sp[j][i]/(float)M);
        
        putchar('\n');
        freq += df;
    }

    for(j = 0; j < 4; j++)
        free(sp[j]);

    return 0;
}
