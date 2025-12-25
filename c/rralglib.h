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

// TERMA DEFINES //
// #define EVENT_BUF_SIZE ((SAMPLE_RATE / 5) * 4)        // sr * 0.8, terma_corrected
// #define CYCLE_BUF_SIZE ((SAMPLE_RATE * 19) + (SAMPLE_RATE/2)) / 5  // sr * 3.9, terma_corrected
#define EVENT_BUF_SIZE (SAMPLE_RATE / 5)        // sr * 0.2, terma_corrected
#define CYCLE_BUF_SIZE ((SAMPLE_RATE * 2) + (SAMPLE_RATE / 10))// sr * 2.1, terma_corrected
// The c preprocessor does not perform floating-point math so this is the workaround to still have these values evaluated at compile time

#define MAIN_BUF_SIZE (SAMPLE_RATE * DATA_WINDOW)
#define MAIN_BUF_HALF_SIZE (SAMPLE_RATE * DATA_WINDOW / 2)
#define MAX_PEAK_COUNT 256

// CWT DEFINES //
#if RRAL_CWT_ENABLED == 1
#include "dsps_fft2r.h"

#define KERNEL_SIZE (SAMPLE_RATE * 3)
#define CWT_SCALES_SIZE 5

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
#if KERNEL_SIZE * 2 - 1 <= 256
#define CWT_OA_SIZE 512
#elif KERNEL_SIZE * 2 - 1 <= 512
#define CWT_OA_SIZE 1024
#elif KERNEL_SIZE * 2 - 1 <= 1024
#define CWT_OA_SIZE 2048
#elif KERNEL_SIZE * 2 - 1 <= 2048
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

#define RRAL_BIQUAD 6
#define SOS_SECT 2
#define SOS_SIZE (SOS_SECT * RRAL_BIQUAD)

#define RUNNING_MEAN_N DATA_WINDOW * SAMPLE_RATE * 3

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

// ALGORITHM STATE STRUCTS
typedef struct state_srmac_t{
    int64_t prevFast;
    int64_t prevSlow;
    int64_t prevCross;
} State_srmac_t;

// JOINT ALGORITHM RESULT STRUCT
typedef struct results_joint_t{
    int32_t hr;
    int32_t rr;
    int32_t sqihr;
    int32_t sqirr;
    int32_t time;
} Results_joint_t;

// JOINT ALGORITHM PARAMETER STRUCT
typedef struct params_joint_t{
    Params_srmac_t srmacParamsHR;
    Params_srmac_t srmacParamsRR;
    int32_t srmacWidthHR;
    int32_t srmacWidthRR;
    int32_t lpfsosHR[SOS_SIZE];
    int32_t lpfsosRR[SOS_SIZE];
    int32_t hpfsosHR[SOS_SIZE];
    int32_t hpfsosRR[SOS_SIZE];
    int32_t runningMeanConstN;
} Params_joint_t;

// JOINT ALGORITHM STATE STRUCT
typedef struct state_joint_t{
    int32_t dataRR[MAIN_BUF_SIZE];
    int32_t dataHR[MAIN_BUF_SIZE];
    uint16_t peaks[MAX_PEAK_COUNT];
    State_srmac_t srmacStateHR;
    State_srmac_t srmacStateRR;
    int32_t hpfxHR[3*SOS_SECT];
    int32_t hpfyHR[3*SOS_SECT];
    int32_t lpfxHR[3*SOS_SECT];
    int32_t lpfyHR[3*SOS_SECT];
    int32_t hpfxRR[3*SOS_SECT];
    int32_t hpfyRR[3*SOS_SECT];
    int32_t lpfxRR[3*SOS_SECT];
    int32_t lpfyRR[3*SOS_SECT];
    int32_t deltaSampleCount;
    int32_t totalSampleCount;
    int64_t runningMeanSum;
    int16_t runningMean;
    int32_t runningMeanN;
} State_joint_t;

// FIXED POINT OPERATIONS
int32_t rral_int_to_fixed(int32_t a);

int32_t rral_fixed_to_int(int32_t a);

float rral_fixed_to_float(int32_t a);

int32_t rral_float_to_fixed(float a);

int32_t rral_mul32(int32_t a, int32_t b);

int32_t rral_div32(int32_t a, int32_t b);

uint32_t rral_isqrt32(uint32_t a);

// FIXED POINT MACROS
#define rral_m_int_to_fixed( Number, Offset )  ((Number) << (Offset))

#define rral_m_fixed_to_int( Number, Offset )  ((Number) >> (Offset))

#define rral_m_fixed_to_float( Number, Offset )  ((float)(Number) / (float)(1 << (Offset)))

#define rral_m_float_to_fixed( Number, Offset )  ((int32_t)((Number) * (float)(1 << (Offset)) + ((Number) >= 0 ? 0.5 : -0.5)))

// FILTERS
int32_t rral_sosfilt(int32_t* data, uint16_t dataLen, int32_t* sos, int32_t sections, int32_t* x, int32_t* y);

// DATA NORMALIZATION FUNCTIONS
void rral_znorm_int(int32_t* data, uint16_t dataLen);

int32_t rral_znorm_fixed(int32_t* data, uint16_t dataLen);

int32_t rral_mean_and_sdev(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* stdev);

int32_t rral_mean_and_sdev_float(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* sdev);

// HELPER FUNCTIONS
uint16_t rral_local_maxima_fast(uint16_t* peaks, uint16_t peaksLen, int32_t* data, uint16_t dataLen);

int32_t rral_arr_min(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound);

int32_t rral_arr_max(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound);

uint8_t rral_remove_by_prominence(uint16_t* peaks, uint16_t peaksLen, uint16_t peakCount, int32_t* data, uint16_t dataLen, int32_t pkProm, int32_t pkHeval, int32_t pkWidth);

uint8_t rral_remove_close_peaks(uint16_t* peaks, uint16_t peaksLen, int16_t peakCount, int32_t* data, uint16_t dataLen, int32_t dist);

// FIND PEAKS ALGORITHM
uint16_t rral_find_peaks(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, Params_findpeaks_t* params);

// TERMA ALGORITHM
int32_t rral_terma(int32_t* data, uint16_t dataLen, int32_t bCoef);

// SRMAC ALGORITHM
int32_t rral_srmac(int32_t* data, uint16_t dataLen, Params_srmac_t params, State_srmac_t *state);

// ZERO-CROSSING BREATH FINDING METHOD
int32_t rral_zero_crossing(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th);

int32_t rral_zero_crossing_raw(int32_t* data, int32_t* rawData, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th);

// RESPIRATORY RATE CALCULATION FUNCTIONS
uint16_t rral_get_respiration_window(uint16_t* peaks, uint16_t peaksLen);

int32_t rral_get_rr(uint16_t* peaks, uint16_t peakCount);

int32_t rral_get_dist(uint16_t* peaks, uint16_t peakCount);

int32_t rral_get_rr_dist(uint16_t* peaks, uint16_t peakCount);

int32_t rral_get_sqi_lite(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount, int32_t mindist, int32_t maxdist);

int32_t rral_get_sqi_full(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount);

int32_t rral_get_sqi_corr(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount);

int32_t rral_joint_hr_rr(int16_t sample, Params_joint_t* params, State_joint_t* state, Results_joint_t* out);

#if RRAL_CWT_ENABLED == 1
float rral_mexh(float t);

float rral_gaus2(float t);

float rral_mexh_scaled(float t, float scale);

float rral_gaus2_scaled(float t, float scale);

int32_t rral_cwt(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table);

int32_t rral_cwt_oa(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table);

int32_t rral_cwt_fixed(int16_t* data, uint16_t dataLen, int16_t minFreq, int16_t maxFreq, int16_t* fft_table);
#endif

#ifdef __cplusplus
}
#endif

#endif /* INC_RRALGLIB_H_ */
