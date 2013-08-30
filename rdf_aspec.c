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
const char version[] = "1.0";

static float lut1bit[256][8];
static float lut2bit[256][4];

static void initluts()
{
    const float HiMag = 3.3359;
	unsigned b, i, l, s, m;
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
 */
static int decode_1bit(const uint8_t * const in, float **out, size_t n)
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

    return 0;
}

/**
 *  parse_rdf_header 
 *  Parses RDF-file header
 *  Print experiment name and date to stdout 
 *  Checks the header size
 *
 *  Returns the number of bytes to offset from the file start to read data
 *  offset = header_size + time_offset
 */
static off_t parse_rdf_header(const char *h, double t_off)
{
    off_t file_header_size = 0;
    char experiment_name[11], data_date[18];
    struct tm tm0;
    time_t t;
    char *time_str;

    strncpy(data_date, &h[6], 18);        
    data_date[17] = 0;

    strncpy(experiment_name, &h[46], 11); 
    experiment_name[10] = 0;

    memset(&tm0, 0, sizeof(tm0));
    strptime(data_date, "%Y %j-%H:%M:%S", &tm0);
    t = mktime(&tm0) + (time_t)t_off;

    file_header_size = (off_t)((unsigned)(h[5] << 8) | h[4]);
    if(file_header_size != 256){
        fprintf(stderr, "file_header_size != 256. Exiting\n");
        exit(EXIT_FAILURE);
    }

    time_str = ctime(&t);
    time_str[strlen(time_str)-1] = 0;
    printf("#%s (%sUT)\n", experiment_name, time_str);

    return file_header_size + (off_t)(400. * 40000. * t_off);
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
    char header[256];
    
    if(argc < 4){
        usage();

        return EXIT_FAILURE;
    }

    f = fopen(argv[1], "r");
    if(f == NULL){
        fprintf(stderr, "Could not open file '%s': ", argv[1]);
        handle_error("");
    }

    /* Number of spectral channels */
    N = atol(argv[2]);
    /* Number of spectra to summ */
    M = atol(argv[3]);

    if(argc == 5)
        time_off = atof(argv[4]);

    n = fread(header, sizeof(uint8_t), 256UL, f);
    if(n != 256){
        fprintf(stderr, "Read only %lu of %lu bytes from RDF-file header. Exiting.", 
                n, 256UL);
        fclose(f);
        exit(EXIT_FAILURE);
    }

    offset = parse_rdf_header(header, time_off);

    if(fseeko(f, offset, SEEK_SET) < 0)
        handle_error("fseeko");

    raw_data = (uint8_t *)malloc(sizeof(uint8_t)*N/2);

    for(j = 0; j < 4; j++){
        sp[j] = calloc(N/2+1, sizeof(float));
        /*c_data[j] = fftwf_alloc_complex(N/2+1);*/ /* Available only since FFTW 3.3-beta1 */
        c_data[j] = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (N/2+1));
        f_data[j] = (float *)c_data[j];
    }

    if(fftwf_import_system_wisdom() != 1)
        fprintf(stderr, "Warning: could not load system fftw wisdom\n");

    p = fftwf_plan_dft_r2c_1d(N, f_data[0], c_data[0], FFTW_ESTIMATE);

    initluts();

    for(k = 0; k < M; k++){
        n = fread(raw_data, sizeof(uint8_t), N/2, f);
        if(n != N/2){
            perror("fread");
            break;
        }
        
        decode_1bit(raw_data, f_data, N);

        /* Four channels */
        for(j = 0; j < 4; j++){
            fftwf_execute_dft_r2c(p, f_data[j], c_data[j]);
            for(i = 0; i <= N/2; i++){
                re = creal(c_data[j][i]);
                im = cimag(c_data[j][i]);
                sp[j][i] += (re*re + im*im)/(float)(N);
            }
        }
    }
    
    free(raw_data);
    fclose(f);

    /* Clean fftw structures */
    fftwf_destroy_plan(p);
    for(j = 0; j < 4; j++){
        fftwf_free(c_data[j]);
    }
    fftwf_cleanup();

    freq = 0.;
    df = 1. / (1e6 * 31.25e-9 * (double)N);

    for(i = 0; i <= N/2; i++){
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
