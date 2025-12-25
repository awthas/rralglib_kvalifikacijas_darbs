#ifndef INC_RRALGLIB_CONFIG_H_
#define INC_RRALGLIB_CONFIG_H_

#ifdef __cplusplus
 extern "C" {
#endif

// SAMPLE RATE IN HZ // Default: 25
#define SAMPLE_RATE 64

// DATA WINDOW SIZE IN SECONDS // Default: 20
#define DATA_WINDOW 20

// THE AMOUNT OF SAMPLES TO ACCUMULATE BEFORE RERUNNING ALGORITHMS // Default: SAMPLE_RATE
#define RRAL_DELTA_SAMPLES SAMPLE_RATE

// LIBRARY FIXED POINT VALUE // Default: 14
#define FIXED_POINT 16

// FILTER COEFFICIENT FIXED POINT VALUE // Default: 20
#define FILT_FIXED_POINT 20

// ENABLE OR DISABLE CWT FUNCTIONS SO THE PROJECT CAN BE COMPILED WITHOUT dsps_fft2r.h // Default: 0
#define RRAL_CWT_ENABLED 0

// MAXIMUM BREATH WIDTH IN SECONDS // Default: 3 // A lower value can lead to less accurate results but makes them more sensitive to lack of detected breathing
#define RRAL_BREATH_WIDTH 3

// MINIMUM RESPIRATORY RATE IN BREATHS PER MINUTE // Default: 6
#define RRAL_MIN_RR 3

#ifdef __cplusplus
}
#endif

#endif