#include "main.h"
#include "unity.h"

void setUp(void) {

}

void tearDown(void) {

}

// rral_int_to_fixed //

void test_ntq_0(void) {
    uint32_t ret = rral_int_to_fixed(0x00000000);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_ntq_1(void) {
    uint32_t ret = rral_int_to_fixed(0x00000001);
    TEST_ASSERT_EQUAL(0x00010000, ret);
}

void test_ntq_65535(void) {
    uint32_t ret = rral_int_to_fixed(0x0000FFFF);
    TEST_ASSERT_EQUAL(0xFFFF0000, ret);
}

void test_ntq_65536(void) {
    uint32_t ret = rral_int_to_fixed(0x00010000);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

// rral_fixed_to_int //

void test_qtn_0_0(void) {
    uint32_t ret = rral_fixed_to_int(0x00000000);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_qtn_1_0(void) {
    uint32_t ret = rral_fixed_to_int(0x00010000);
    TEST_ASSERT_EQUAL(0x00000001, ret);
}

void test_qtn_1_1(void) {
    uint32_t ret = rral_fixed_to_int(0x00018000);
    TEST_ASSERT_EQUAL(0x00000001, ret);
}

void test_qtn_65535_65535(void) {
    uint32_t ret = rral_fixed_to_int(0xFFFFFFFF);
    TEST_ASSERT_EQUAL(0xFFFFFFFF, ret);
}

// rral_float_to_fixed //

void test_ftq_0_0(void) {
    uint32_t ret = rral_float_to_fixed(0.0f);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_ftq_1_1(void) {
    uint32_t ret = rral_float_to_fixed(1.5f);
    TEST_ASSERT_EQUAL(0x00018000, ret);
}

void test_ftq_32767_0(void) {
    uint32_t ret = rral_float_to_fixed(32767.0f);
    TEST_ASSERT_EQUAL(0x7FFF0000, ret);
}

void test_ftq_32768_0(void) {
    uint32_t ret = rral_float_to_fixed(32768.0f);
    TEST_ASSERT_EQUAL(0x80000000, ret);
}

// rral_fixed_to_float //

void test_qtf_0_0(void) {
    float ret = rral_fixed_to_float(0x00000000);
    TEST_ASSERT_EQUAL_FLOAT(0.0f, ret);
}

void test_qtf_1_1(void) {
    float ret = rral_fixed_to_float(0x00018000);
    TEST_ASSERT_EQUAL_FLOAT(1.5f, ret);
}

void test_qtf_0_65534(void) {
    float ret = rral_fixed_to_float(0x0000FFFE);
    TEST_ASSERT_EQUAL_FLOAT(0.999969f, ret);
}

void test_qtf_32767_65535(void) {
    float ret = rral_fixed_to_float(0x7FFFFFFF);
    TEST_ASSERT_EQUAL_FLOAT(32768.0f, ret);
}

void test_qtf_32767_0(void) {
    float ret = rral_fixed_to_float(0x7FFF0000);
    TEST_ASSERT_EQUAL_FLOAT(32767.0f, ret);
}

// rral_mul32 //

void test_mul32_zero(void) {
    int32_t ret = rral_mul32(0, 1);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_mul32_regular(void) {
    int32_t ret = rral_mul32(32767, 32767);
    TEST_ASSERT_EQUAL((0x3FFF0001 >> FIXED_POINT), ret);
}

void test_mul32_ovf(void) {
    int32_t ret = rral_mul32(65536, 65536);
    TEST_ASSERT_EQUAL(0x100000000 >> FIXED_POINT, ret);
}

// rral_div32 //

void test_div32_zero(void) {
    int32_t ret = rral_div32(1, 0);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_div32_one(void) {
    int32_t ret = rral_div32(0x00010000, 0x00000001);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_div32_one_q(void) {
    int32_t ret = rral_div32(0x00010000, 0x00010000);
    TEST_ASSERT_EQUAL(0x00010000, ret);
}

// rral_isqrt32 // 

void test_isqrt32_zero(void) {
    int32_t ret = rral_isqrt32(0);
    TEST_ASSERT_EQUAL(0x00000000, ret);
}

void test_isqrt32_one(void) {
    int32_t ret = rral_isqrt32(1);
    TEST_ASSERT_EQUAL(0x00000001, ret);
}

void test_isqrt32_max(void) {
    int32_t ret = rral_isqrt32(0xFFFFFFFF);
    TEST_ASSERT_EQUAL(0x0000FFFF, ret); // 65535.999... rounded down
}

// rral_mean_and_stdev //

void test_meanstdev_zeroes(void) {
    int32_t data[16] = { 0 };
    int32_t mean = 0;
    int32_t stdev = 0;

    int32_t ret = rral_mean_and_sdev(data, 16, &mean, &stdev);

    TEST_ASSERT_EQUAL(0x00000000, mean);
    TEST_ASSERT_EQUAL(0x00000000, stdev);
}

void test_meanstdev_ones(void) {
    int32_t data[16] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    int32_t mean = 0;
    int32_t stdev = 0;

    int32_t ret = rral_mean_and_sdev(data, 16, &mean, &stdev);

    TEST_ASSERT_EQUAL(0x00000001, mean);
    TEST_ASSERT_EQUAL(0x00000000, stdev);
}

void test_meanstdev_sine(void) {
    // sin[0:2pi]
    int32_t data[16] = { 0, 406, 743, 951, 994, 866, 587, 207, -207, -587, -866, -994, -951, -743, -406, 0 };
    int32_t mean = 0;
    int32_t stdev = 0;

    for (uint16_t i = 0; i < 16; i++) {
        data[i] = rral_int_to_fixed(data[i]);
    }

    int32_t ret = rral_mean_and_sdev(data, 16, &mean, &stdev);

    TEST_ASSERT_EQUAL(0x00000000, mean);
    TEST_ASSERT_EQUAL(0x02C20000, (stdev &= 0xFFFF0000)); // A fair amount of error is permissible so check only if the integer part matches
    //TEST_ASSERT_EQUAL_FLOAT(706.7f, rral_fixed_to_float(stdev));
}

void test_meanstdev_unitvar(void) {
    // sin[0:2pi]
    int32_t data[16] = {
        rral_float_to_fixed(0.0f),
        rral_float_to_fixed(0.5744238149250852f),
        rral_float_to_fixed(1.0512238780525573f),
        rral_float_to_fixed(1.3455099704279705f),
        rral_float_to_fixed(1.4063479606786569f),
        rral_float_to_fixed(1.2252488269091717f),
        rral_float_to_fixed(0.8305093087709975f),
        rral_float_to_fixed(0.292871255392839f),
        rral_float_to_fixed(-0.292871255392839f),
        rral_float_to_fixed(-0.8305093087709975f),
        rral_float_to_fixed(-1.2252488269091717f),
        rral_float_to_fixed(-1.4063479606786569f),
        rral_float_to_fixed(-1.3455099704279705f),
        rral_float_to_fixed(-1.0512238780525573f),
        rral_float_to_fixed(-0.5744238149250852f),
        rral_float_to_fixed(0.0f)
    };
    int32_t mean = 0;
    int32_t stdev = 0;

    int32_t ret = rral_mean_and_sdev(data, 16, &mean, &stdev);

    TEST_ASSERT_EQUAL(0x00000000, mean);
    TEST_ASSERT_LESS_OR_EQUAL_FLOAT(1.01f, rral_fixed_to_float(stdev)); // A fair amount of error is permissible
    TEST_ASSERT_GREATER_OR_EQUAL_FLOAT(0.99f, rral_fixed_to_float(stdev));
}

void test_meanstdev_noise(void) {
    int32_t data[16] = {
        rral_int_to_fixed(15041),
        rral_int_to_fixed(7041),
        rral_int_to_fixed(-31479),
        rral_int_to_fixed(-14185),
        rral_int_to_fixed(12131),
        rral_int_to_fixed(-24782),
        rral_int_to_fixed(3596),
        rral_int_to_fixed(5333),
        rral_int_to_fixed(-27129),
        rral_int_to_fixed(-23293),
        rral_int_to_fixed(30407),
        rral_int_to_fixed(16953),
        rral_int_to_fixed(1779),
        rral_int_to_fixed(22948),
        rral_int_to_fixed(13015),
        rral_int_to_fixed(-12641)
    };

    int32_t mean = 0;
    int32_t stdev = 0;

    int32_t ret = rral_mean_and_sdev(data, 16, &mean, &stdev);

    TEST_ASSERT_EQUAL_FLOAT(-329.0625f, rral_fixed_to_float(mean));
    TEST_ASSERT_EQUAL_FLOAT(19353.54f, rral_fixed_to_float(stdev));
}

void test_sqifull_pattern(void) {
    int32_t pattern[8] = {
        rral_int_to_fixed(1),
        rral_int_to_fixed(2),
        rral_int_to_fixed(3),
        rral_int_to_fixed(4),
        rral_int_to_fixed(5),
        rral_int_to_fixed(4),
        rral_int_to_fixed(3),
        rral_int_to_fixed(2),
    };
    int32_t data[128];

    for (int32_t i = 0; i < 128; i++) {
        data[i] = pattern[i % 8];
    }

    uint16_t peaks[4] = { 20, 28, 36, 52 };

    int32_t sqi = rral_get_sqi_full(data, 128, peaks, 4);
    TEST_ASSERT_EQUAL_FLOAT(1.0f, rral_fixed_to_float(sqi));
}

void test_sqifull_twopeaks(void) {
    int32_t pattern[8] = {
        rral_int_to_fixed(1),
        rral_int_to_fixed(2),
        rral_int_to_fixed(3),
        rral_int_to_fixed(4),
        rral_int_to_fixed(5),
        rral_int_to_fixed(4),
        rral_int_to_fixed(3),
        rral_int_to_fixed(2),
    };
    int32_t data[128];

    for (int32_t i = 0; i < 128; i++) {
        data[i] = pattern[i % 8];
    }

    uint16_t peaks[2] = { 24, 48 };

    int32_t sqi = rral_get_sqi_full(data, 128, peaks, 2);
    TEST_ASSERT_EQUAL_FLOAT(1.0f, rral_fixed_to_float(sqi));
}

void test_arrmin_one(void) {
    int32_t arr[1] = { 0 };
    TEST_ASSERT_EQUAL(0, rral_arr_min(arr, 1, 0, 1));
}

void test_arrmin_zeroes(void) {
    int32_t arr[8] = { 0,0,0,0,0,0,0,0 };
    TEST_ASSERT_EQUAL(0, rral_arr_min(arr, 8, 0, 8));
}

void test_arrmin_twotroughs(void) {
    int32_t arr[16] = { 0,0,0,0,-1,0,0,0,0,0,0,0,0,-1,0,0 };
    TEST_ASSERT_EQUAL(-1, rral_arr_min(arr, 16, 0, 16));
}

void test_arrmax_one(void) {
    int32_t arr[1] = { 0 };
    TEST_ASSERT_EQUAL(0, rral_arr_max(arr, 1, 0, 1));
}

void test_arrmax_zeroes(void) {
    int32_t arr[8] = { 0,0,0,0,0,0,0,0 };
    TEST_ASSERT_EQUAL(0, rral_arr_max(arr, 8, 0, 8));
}

void test_arrmax_twopeaks(void) {
    int32_t arr[16] = { 0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0 };
    TEST_ASSERT_EQUAL(1, rral_arr_max(arr, 16, 0, 16));
}

void test_zerocross_zeroes(void) {
    int32_t data[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    uint16_t peaks[2] = { 0, 0 };

    int32_t count = rral_zero_crossing(data, 16, peaks, 2, 1, 0);

    TEST_ASSERT_EQUAL(0, count);
    TEST_ASSERT_EQUAL(0, peaks[0]);
}

void test_zerocross_one(void) {
    int32_t data[1] = { 1 };
    uint16_t peaks[2] = { 0, 0 };

    int32_t count = rral_zero_crossing(data, 1, peaks, 2, 1, 0);

    TEST_ASSERT_EQUAL(0, count);
    TEST_ASSERT_EQUAL(0, peaks[0]);
}

void test_zerocross_impulse(void) {
    int32_t data[16] = { 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0 };
    uint16_t peaks[2] = { 0, 0 };

    int32_t count = rral_zero_crossing(data, 16, peaks, 2, 1, 0);

    TEST_ASSERT_EQUAL(1, count);
    TEST_ASSERT_EQUAL(5, peaks[0]);
}

void test_zerocross_sine(void) {
    // sin[0:2pi]
    int32_t data[16] = { 0, 406, 743, 951, 994, 866, 587, 207, -207, -587, -866, -994, -951, -743, -406, 0 };
    uint16_t peaks[2] = { 0, 0 };

    int32_t count = rral_zero_crossing(data, 16, peaks, 2, 1, 0);

    TEST_ASSERT_EQUAL(1, count);
    TEST_ASSERT_EQUAL(4, peaks[0]);
}

void test_zerocrossraw_sineswap(void) {
    // sin[0:2pi] // raw data has peak in a slightly different spot
    int32_t data[16] = { 0, 406, 743, 951, 994, 866, 587, 207, -207, -587, -866, -994, -951, -743, -406, 0 };
    int32_t dataraw[16] = { 0, 406, 743, 994, 951, 866, 587, 207, -207, -587, -866, -994, -951, -743, -406, 0 };
    uint16_t peaks[2] = { 0, 0 };

    int32_t count = rral_zero_crossing_raw(data, dataraw, 16, peaks, 2, 1, 0);

    TEST_ASSERT_EQUAL(1, count);
    TEST_ASSERT_EQUAL(3, peaks[0]);
}

void test_findpeaks_zeroes(void) {
    int32_t data[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    uint16_t peaks[2] = { 0, 0 };
    Params_findpeaks_t params = {1, 0, 0, 0, 0};

    int32_t count = rral_find_peaks(data, 16, peaks, 2, &params);

    TEST_ASSERT_EQUAL(0, count);
    TEST_ASSERT_EQUAL(0, peaks[0]);
}

void test_findpeaks_one(void) {
    int32_t data[1] = { 1 };
    uint16_t peaks[2] = { 0, 0 };
    Params_findpeaks_t params = { 1, 0, 0, 0, 0 };

    int32_t count = rral_find_peaks(data, 1, peaks, 2, &params);

    TEST_ASSERT_EQUAL(0, count);
    TEST_ASSERT_EQUAL(0, peaks[0]);
}

void test_findpeaks_impulse(void) {
    int32_t data[16] = { 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0 };
    uint16_t peaks[2] = { 0, 0 };
    Params_findpeaks_t params = { 1, 0, 0, 0, 0 };

    int32_t count = rral_find_peaks(data, 16, peaks, 2, &params);

    TEST_ASSERT_EQUAL(1, count);
    TEST_ASSERT_EQUAL(5, peaks[0]);
}

void test_findpeaks_sine(void) {
    // sin[0:2pi]
    int32_t data[16] = { 0, 406, 743, 951, 994, 866, 587, 207, -207, -587, -866, -994, -951, -743, -406, 0 };
    uint16_t peaks[2] = { 0, 0 };
    Params_findpeaks_t params = { 1, 0, 0, 0, 0 };

    int32_t count = rral_find_peaks(data, 16, peaks, 2, &params);

    TEST_ASSERT_EQUAL(1, count);
    TEST_ASSERT_EQUAL(4, peaks[0]);
}

int main()
{
    char temp[32];
    uint16_t templ = 0;
    templ = sprintf_s(temp, 32, "--- RRALGLIB UNIT TESTING ---\n");
    printf(temp);

    UNITY_BEGIN();

    // rral_int_to_fixed
    templ = sprintf_s(temp, 32, "\n- rral_int_to_fixed -\n");
    printf(temp);
    RUN_TEST(test_ntq_0);
    RUN_TEST(test_ntq_1);
    RUN_TEST(test_ntq_65535);
    RUN_TEST(test_ntq_65536);

    // rral_fixed_to_int
    templ = sprintf_s(temp, 32, "\n- rral_fixed_to_int -\n");
    printf(temp);
    RUN_TEST(test_qtn_0_0);
    RUN_TEST(test_qtn_1_0);
    RUN_TEST(test_qtn_1_1);
    RUN_TEST(test_qtn_65535_65535);

    // rral_float_to_fixed
    templ = sprintf_s(temp, 32, "\n- rral_float_to_fixed -\n");
    printf(temp);
    RUN_TEST(test_ftq_0_0);
    RUN_TEST(test_ftq_1_1);
    RUN_TEST(test_ftq_32767_0);
    RUN_TEST(test_ftq_32768_0);
    
    // rral_fixed_to_float
    templ = sprintf_s(temp, 32, "\n- rral_fixed_to_float -\n");
    printf(temp);
    RUN_TEST(test_qtf_0_0);
    RUN_TEST(test_qtf_1_1);
    RUN_TEST(test_qtf_0_65534);
    RUN_TEST(test_qtf_32767_65535);
    RUN_TEST(test_qtf_32767_0);

    // rral_mul32
    templ = sprintf_s(temp, 32, "\n- rral_mul32 -\n");
    printf(temp);
    RUN_TEST(test_mul32_zero);
    RUN_TEST(test_mul32_regular);
    RUN_TEST(test_mul32_ovf);

    // rral_div32
    templ = sprintf_s(temp, 32, "\n- rral_div32 -\n");
    printf(temp);
    RUN_TEST(test_div32_zero);
    RUN_TEST(test_div32_one);
    RUN_TEST(test_div32_one_q);

    // rral_isqrt32
    templ = sprintf_s(temp, 32, "\n- rral_isqrt32 -\n");
    printf(temp);
    RUN_TEST(test_isqrt32_zero);
    RUN_TEST(test_isqrt32_one);
    RUN_TEST(test_isqrt32_max);
    
    // rral_mean_and_stdev
    templ = sprintf_s(temp, 32, "\n- rral_mean_and_stdev -\n");
    printf(temp);
    RUN_TEST(test_meanstdev_zeroes);
    RUN_TEST(test_meanstdev_ones);
    RUN_TEST(test_meanstdev_sine);
    RUN_TEST(test_meanstdev_unitvar);
    RUN_TEST(test_meanstdev_noise);

    // rral_arrmin
    templ = sprintf_s(temp, 32, "\n- rral_arr_min -\n");
    printf(temp);
    RUN_TEST(test_arrmin_one);
    RUN_TEST(test_arrmin_zeroes);
    RUN_TEST(test_arrmin_twotroughs);

    // rral_arrmax
    templ = sprintf_s(temp, 32, "\n- rral_arr_max -\n");
    printf(temp);
    RUN_TEST(test_arrmax_one);
    RUN_TEST(test_arrmax_zeroes);
    RUN_TEST(test_arrmax_twopeaks);

    // rral_zero_crossing
    templ = sprintf_s(temp, 32, "\n- rral_zero_crossing -\n");
    printf(temp);
    RUN_TEST(test_zerocross_zeroes);
    RUN_TEST(test_zerocross_one);
    RUN_TEST(test_zerocross_impulse);
    RUN_TEST(test_zerocross_sine);

    // rral_zero_crossing_raw
    templ = sprintf_s(temp, 32, "\n- rral_zero_crossing_raw -\n");
    printf(temp);
    RUN_TEST(test_zerocrossraw_sineswap);

    // rral_find_peaks
    templ = sprintf_s(temp, 32, "\n- rral_find_peaks -\n");
    printf(temp);
    RUN_TEST(test_findpeaks_zeroes);
    RUN_TEST(test_findpeaks_one);
    RUN_TEST(test_findpeaks_impulse);
    RUN_TEST(test_findpeaks_sine);

    // rral_srmac
    //templ = sprintf_s(temp, 32, "\n- rral_srmac -\n");
    //printf(temp);
    //RUN_TEST(test_srmac_zeroes);
    //RUN_TEST(test_srmac_one);
    //RUN_TEST(test_srmac_impulse);
    //RUN_TEST(test_srmac_sine);

    // rral_terma
    //templ = sprintf_s(temp, 32, "\n- rral_terma -\n");
    //printf(temp);
    //RUN_TEST(test_terma_zeroes);
    //RUN_TEST(test_terma_one);
    //RUN_TEST(test_terma_impulse);
    //RUN_TEST(test_terma_sine);

    // rral_wavelet



    // rral_cwt_oa



    // rral_sqi_full
    templ = sprintf_s(temp, 32, "\n- rral_sqi_full -\n");
    printf(temp);
    RUN_TEST(test_sqifull_pattern);
    RUN_TEST(test_sqifull_twopeaks);

	return UNITY_END();
}
