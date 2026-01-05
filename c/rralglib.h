#ifndef INC_RRALGLIB_H_
#define INC_RRALGLIB_H_

#ifdef __cplusplus
 extern "C" {
#endif

// LIBRARY CONFIGURATION FILE //
#include "rralglib_config.h"

// STANDARD INCLUDES //
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

// MAIN BUFFER SIZE //
#if (SAMPLE_RATE * DATA_WINDOW) % 2 == 0
#define MAIN_BUF_SIZE (SAMPLE_RATE * DATA_WINDOW)
#else
#define MAIN_BUF_SIZE ((SAMPLE_RATE * DATA_WINDOW) + 1)
#endif

// HALF OF THE MAIN BUFFER //
#define MAIN_BUF_HALF_SIZE (MAIN_BUF_SIZE / 2)

// MAXIMUM PEAK COUNT //
#define MAX_PEAK_COUNT 256

// TERMA DEFINES //
#define EVENT_BUF_SIZE (SAMPLE_RATE / 10)
#define CYCLE_BUF_SIZE (SAMPLE_RATE / 20)

// CWT DEFINES //
#if RRAL_CWT_ENABLED == 1
#include "dsps_fft2r.h"

// KERNEL SIZE IN SAMPLES //
#define KERNEL_SIZE (SAMPLE_RATE * 3)

// CWT_BUF_SIZE DEFINITION BASED ON THE INPUT DATA SIZE AND THE KERNEL SIZE //
#ifndef CWT_BUF_SIZE
#if ((SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) <= 256)
#define CWT_BUF_SIZE 256
#elif ((SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) <= 512)
#define CWT_BUF_SIZE 512
#elif ((SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) <= 1024)
#define CWT_BUF_SIZE 1024
#elif ((SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) <= 2048)
#define CWT_BUF_SIZE 2048
#elif ((SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) <= 4096)
#define CWT_BUF_SIZE 4096
#else /* (SAMPLE_RATE * DATA_WINDOW + KERNEL_SIZE - 1) > 2048 */
#define CWT_BUF_SIZE 8192
#endif
#endif

// CWT_OA_SIZE DEFINITION
#ifndef CWT_OA_SIZE
#if (KERNEL_SIZE * 2 - 1) <= 256
#define CWT_OA_SIZE 512
#elif (KERNEL_SIZE * 2 - 1) <= 512
#define CWT_OA_SIZE 1024
#elif (KERNEL_SIZE * 2 - 1) <= 1024
#define CWT_OA_SIZE 2048
#elif (KERNEL_SIZE * 2 - 1) <= 2048
#define CWT_OA_SIZE 4096
#else
#define CWT_OA_SIZE 8192
#endif
#endif

// DEFAULT CWT_BUF_SIZE FALLBACK VALUE //
#ifndef CWT_BUF_SIZE
#define CWT_BUF_SIZE 512*2
#endif
#endif

// BIQUAD SIZE //
#define RRAL_BIQUAD 6

// BIQUAD SECTION COUNT (FILTER ORDER - 1) //
#define SOS_SECT 2

// FULL SOS COEFFICIENT ARRAY SIZE
#define SOS_SIZE (SOS_SECT * RRAL_BIQUAD)

// REALTIME RUNNING MEAN MAXIMUM WINDOW SIZE IN SAMPLES //
#define RUNNING_MEAN_N (DATA_WINDOW * SAMPLE_RATE * 3)

// RET ENUM
typedef enum {
	RRAL_OK = 0,
	RRAL_ERR = -1,
	RRAL_RDY = 1
} RRAL_Status_t;

// ALGORITHM PARAMETER STRUCTS
typedef struct params_srmac_t{
    int32_t coefFast;
    int32_t coefSlow;
    int32_t coefCross;
} Params_srmac_t;

typedef struct params_findpeaks_t{
    int32_t pkProm;
    int32_t pkHeval;
    int32_t pkWidth;
    int32_t pkProxim;
    int32_t pkOrder;
} Params_findpeaks_t;

typedef struct params_peaks_t{
    Params_srmac_t srmacParams;
    int32_t srmacWidth;
    int32_t lpfsos[SOS_SIZE];
    int32_t hpfsos[SOS_SIZE];
} Params_peaks_t;

// ALGORITHM STATE STRUCTS
typedef struct state_srmac_t{
    int64_t prevFast;
    int64_t prevSlow;
    int64_t prevCross;
} State_srmac_t;

typedef struct state_peaks_t{
    int32_t data[MAIN_BUF_SIZE];
    uint16_t peaks[MAX_PEAK_COUNT];
    uint16_t peaksN;
    State_srmac_t srmacState;
    int32_t hpfx[3*SOS_SECT];
    int32_t hpfy[3*SOS_SECT];
    int32_t lpfx[3*SOS_SECT];
    int32_t lpfy[3*SOS_SECT];
    int32_t deltaSampleCount;
} State_peaks_t;

typedef struct state_joint_t{
    int32_t deltaSampleCount;
    int32_t totalSampleCount;
    int64_t runningMeanSum;
    int32_t runningMean;
} State_joint_t;

typedef struct results_joint_t{
    int32_t hr;
    int32_t rr;
    int32_t sqihr;
    int32_t sqirr;
    int32_t time;
} Results_joint_t;

typedef struct results_single_t{
    int32_t rate;
    int32_t sqi;
    int32_t time;
} Results_single_t;

// FIXED POINT OPERATIONS
int32_t rral_int_to_fixed(int32_t a);

int32_t rral_fixed_to_int(int32_t a);

float rral_fixed_to_float(int32_t a);

int32_t rral_float_to_fixed(float a);

int32_t rral_mul32(int32_t a, int32_t b);

int32_t rral_div32(int32_t a, int32_t b);

uint32_t rral_isqrt32(uint32_t a);

uint64_t rral_isqrt64(uint64_t a);

// FILTERS
int32_t rral_sosfilt(int32_t* data, uint16_t dataLen, int32_t* sos, int32_t sections, int32_t* x, int32_t* y);

// GENERAL FUNCTIONS
int32_t rral_znorm_fixed(int32_t* data, uint16_t dataLen);

int32_t rral_mean_and_sdev(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* stdev);

int32_t rral_mean_and_sdev_float(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* sdev);

int32_t rral_arr_min(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound);

int32_t rral_arr_max(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound);

// FIND_PEAKS ALGORITHM
int32_t rral_local_maxima_fast(uint16_t* peaks, uint16_t peaksLen, int32_t* data, uint16_t dataLen);

int32_t rral_remove_by_prominence(uint16_t* peaks, uint16_t peaksLen, uint16_t peakCount, int32_t* data, uint16_t dataLen, int32_t pkProm, int32_t pkHeval, int32_t pkWidth);

int32_t rral_remove_close_peaks(uint16_t* peaks, uint16_t peaksLen, int16_t peakCount, int32_t* data, uint16_t dataLen, int32_t dist);

int32_t rral_find_peaks(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, Params_findpeaks_t* params);

// TERMA ALGORITHM
int32_t rral_terma(int32_t* data, uint16_t dataLen, int32_t bCoef);

// SRMAC ALGORITHM
int32_t rral_srmac(int32_t* data, uint16_t dataLen, Params_srmac_t params, State_srmac_t *state);

// ZERO-CROSSING BREATH FINDING METHOD
int32_t rral_zero_crossing(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th);

int32_t rral_zero_crossing_raw(int32_t* data, int32_t* rawData, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th);

// SQI METHODS
int32_t rral_get_sqi_lite(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount, int32_t mindist, int32_t maxdist);

int32_t rral_get_sqi_full(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount);

// REALTIME PROCESSING ALGORITHMS
int32_t rral_state_peaks_init(State_peaks_t* state);

int32_t rral_state_joint_init(State_joint_t* state);

int32_t rral_realtime_peaks(Params_peaks_t* params, State_peaks_t* state);

int32_t rral_realtime_peaks_shift(State_peaks_t* state);

int32_t rral_single_param(int16_t sample, Params_peaks_t* params, State_peaks_t* state, State_joint_t* stateJ, Results_single_t* out);

int32_t rral_joint_hr_rr(int16_t sample, Params_peaks_t* paramsHR, Params_peaks_t* paramsRR, State_peaks_t* stateHR, State_peaks_t* stateRR, State_joint_t* stateJ, Results_joint_t* out);

// CONDITIONAL COMPILATION FOR OPTIONAL FFT-BASED CWT
#if RRAL_CWT_ENABLED == 1
float rral_mexh(float t);

float rral_gaus2(float t);

float rral_mexh_scaled(float t, float scale);

float rral_gaus2_scaled(float t, float scale);

float rral_wavelet_scaled(float t, float scale, uint8_t wavelet);

int32_t rral_cwt(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table);

int32_t rral_cwt_oa(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table);

int32_t rral_cwt_oa_fast_init(float* kernels, int32_t kernelBufLen, float minFreq, float maxFreq, float* fft_table);

int32_t rral_cwt_fixed(int16_t* data, uint16_t dataLen, int16_t minFreq, int16_t maxFreq, int16_t* fft_table);
#endif

#ifdef __cplusplus
}
#endif

#endif /* INC_RRALGLIB_H_ */
