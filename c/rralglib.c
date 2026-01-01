#include "rralglib.h"

// FIXED POINT OPERATIONS
int32_t rral_int_to_fixed(int32_t a){
	return a << FIXED_POINT;
}

int32_t rral_fixed_to_int(int32_t a){
	return (int32_t)(a >> FIXED_POINT);
}

float rral_fixed_to_float(int32_t a){
	return (float)a / (float)(1 << FIXED_POINT);
}

int32_t rral_float_to_fixed(float a){
	return (int32_t)(a * (float)(1 << FIXED_POINT) + (a >= 0 ? 0.5 : -0.5));
}

int32_t rral_mul32(int32_t a, int32_t b){
	return ((int64_t)a*(int64_t)b) >> FIXED_POINT;
}

// UNFINISHED MUL64 FUNCTION
// uint64_t rral_mul64(uint64_t a, uint64_t b){ //https://stackoverflow.com/a/58381061
// 	uint64_t lo_lo = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF);
// 	uint64_t hi_lo = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF);
// 	uint64_t lo_hi = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF);
// 	uint64_t hi_hi = (a & 0xFFFFFFFF) * (b & 0xFFFFFFFF);
// 	return ((int64_t)a*(int64_t)b) >> FIXED_POINT;
// }

int32_t rral_div32(int32_t a, int32_t b){
	if(b == 0){return 0;}
	return (((int64_t)a) << FIXED_POINT) / (int64_t)b;
}

// CORDIC 32-bit integer square root function // https://www.convict.lu/Jeunes/Math/square_root_CORDIC.htm
uint32_t rral_isqrt32(uint32_t a){
	int64_t res = 0;
	int32_t base = 32768;
	for(uint16_t i = 1; i <= 16; i++){
		res += base;
		if((res * res) > a){
			res -= base;
		}
		base = base >> 1;
	}
	return (uint32_t)res;
}

uint64_t rral_isqrt64(uint64_t a){
	int64_t res = 0;
	int64_t base = 2147483648;
	for(uint16_t i = 1; i <= 32; i++){
		res += base;
		if((res * res) > a){
			res -= base;
		}
		base = base >> 1;
	}
	return (uint64_t)res;
}

int32_t sqrt64f(int64_t a){
	double a_f = (double)a / (double)(1 << FIXED_POINT);
	a_f = sqrt(a_f);
	return (int32_t)((float)a_f * (float)(1 << FIXED_POINT) + (a_f >= 0 ? 0.5 : -0.5));
}

// Biquad IIR filter function(sos coefs should have a[0] normalized to 1.0)
int32_t rral_sosfilt(int32_t* data, uint16_t dataLen, int32_t* sos, int32_t sections, int32_t* x, int32_t* y){
	if(dataLen <= 0){ return RRAL_ERR;}
	if(sections <= 0){ return RRAL_ERR;}
	
	// Buffers
	int32_t b[3];
	int32_t a[3];
	int64_t temp = 0; // used in place of y[0]
	int32_t offset = 0; // used to store section offsets

	// Filter data once for each sos section
	for(uint16_t sec = 0; sec < sections; sec++){
		// Grab coefficients
		offset = sec * 6;
		b[0] = sos[offset];
		b[1] = sos[offset+1];
		b[2] = sos[offset+2];
		a[0] = sos[offset+3];
		a[1] = sos[offset+4];
		a[2] = sos[offset+5];

		offset = sec * 3;
		for(int32_t i = 0; i < dataLen; i++){
			x[offset] = data[i];

			// Apply filter
			temp = ((int64_t)b[0]*x[offset] + (int64_t)b[1]*x[offset+1] + (int64_t)b[2]*x[offset+2] - (int64_t)a[1]*y[offset+1] - (int64_t)a[2]*y[offset+2]) >> FIXED_POINT;
			if(temp > 0x7FFFFFFF){ temp = 0x7FFFFFFF;} 		// Clamp if y[0] overflows 32int_t
			// if(temp < -0x80000000){ temp = -0x80000000;} 	// ...
			data[i] = (int32_t)temp;

			// Update buffers
			x[offset+2] = x[offset+1];
			x[offset+1] = x[offset];
			y[offset+2] = y[offset+1];
			y[offset+1] = (int32_t)temp;
		}
	}

	return RRAL_OK;
}

// Mean and standard deviation calculation function using the two-pass algorithm
int32_t rral_mean_and_sdev(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* stdev){
	if(dataLen <= 1){ return RRAL_ERR; }
	if(mean == NULL){ return RRAL_ERR; }
	if(stdev == NULL){ return RRAL_ERR; }

	// Calculate the mean
	int64_t meanAcc = 0;
	for(int32_t i = 0; i < dataLen; i++){
		meanAcc += data[i];
	}
	*mean = (meanAcc / dataLen);

	// Calculate the sum of differences to the mean squared
	int64_t sumOfSquares = 0;
	int64_t delta = 0;
	for(int32_t i = 0; i < dataLen; i++){
		delta = data[i] - *mean;
		sumOfSquares += ((delta * delta) >> FIXED_POINT);
	}

	// Standard deviation is calculated with ddof = 1
	
	// This approach works but is both imprecise and prone to overflows when input data is large or far from the mean
	//*stdev = (int64_t)((sumOfSquares / (dataLen - 1)) >> 8);  // these two steps lose quite a lot of precision but at least it works
	//*stdev = rral_isqrt32(*stdev) << (FIXED_POINT - (8 / 2)); // add a check for FIXED_POINT size

	// For the sake of precision the final step has been implemented using double precision
	sumOfSquares = sumOfSquares / (dataLen - 1);
	double d_stdev = (double)sumOfSquares / (double)(1 << FIXED_POINT);
	d_stdev = sqrt(d_stdev);
	*stdev = (int32_t)((float)d_stdev * (float)(1 << FIXED_POINT) + (d_stdev >= 0 ? 0.5 : -0.5));
		
	return RRAL_OK;
}

int32_t rral_mean_and_sdev_float(int32_t* data, uint16_t dataLen, int32_t* mean, int32_t* sdev){
	if(dataLen <= 1){ return RRAL_ERR; }
	if(mean == NULL){ return RRAL_ERR; }
	if(sdev == NULL){ return RRAL_ERR; }

	// Calculate the mean
	int64_t meanAcc = 0;
	for(int32_t i = 0; i < dataLen; i++){
		meanAcc += data[i];
	}
	*mean = (meanAcc / dataLen);

	// Calculate the sum of differences to the mean squared
	float sumOfSquares = 0;
	float fpconst = (float)(1 << FIXED_POINT);
	int64_t delta = 0;
	for (int32_t i = 0; i < dataLen; i++) {
		delta = data[i] - *mean;
		sumOfSquares += (float)((delta * delta) >> FIXED_POINT) / fpconst;
	}

	// Calculate the standard deviation
	sumOfSquares = sumOfSquares / (dataLen - 1);
	sumOfSquares = sqrtf(sumOfSquares);
	*sdev = rral_float_to_fixed(sumOfSquares);
		
	return RRAL_OK;
}

// Local maxima detection without an order parameter
int32_t rral_local_maxima_fast(uint16_t* peaks, uint16_t peaksLen, int32_t* data, uint16_t dataLen) {
	if (dataLen == 0) { return RRAL_ERR; }

	int32_t lastY = -2147483648;
	uint16_t ind = 0;
	for (uint16_t i = 1; i < dataLen - 1; i++) {
		// Check for plateau
		if (data[i] == lastY) {
			continue;
		}
		lastY = data[i];

		// If this is a local maxima add it to the peaks buffer
		if (data[i - 1] <= data[i] && data[i + 1] <= data[i]) {
			peaks[ind] = i;
			if (ind + 1 < peaksLen) { ind++; }
		}
	}
	return ind;
}


// Return min array value
int32_t rral_arr_min(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound){
	if(left_bound <= 0){
		left_bound = 0;
	}
	if(dataLen <= 0){
		return 0;
	}
	if(right_bound >= dataLen){
		right_bound = dataLen-1;
	}
	if(left_bound > right_bound){
		left_bound = 0;
		right_bound = dataLen-1;
	}
	if(left_bound == right_bound){
		return data[left_bound];
	}
	int32_t min = data[left_bound];
	for(uint16_t i = left_bound+1; i <= right_bound; i++){
		if(data[i] < min){
			min = data[i];
		}
	}
	return min;
}

// Return max array value
int32_t rral_arr_max(int32_t* data, uint16_t dataLen, uint16_t left_bound, uint16_t right_bound){
	if(left_bound <= 0){
		left_bound = 0;
	}
	if(dataLen <= 0){
		return 0;
	}
	if(right_bound >= dataLen){
		right_bound = dataLen-1;
	}
	if(left_bound > right_bound){
		left_bound = 0;
		right_bound = dataLen-1;
	}
	if(left_bound == right_bound){
		return data[left_bound];
	}
	int32_t max = data[left_bound];
	for(uint16_t i = left_bound+1; i <= right_bound; i++){
		if(data[i] > max){
			max = data[i];
		}
	}
	return max;
}

// The optimized version of peak prominence filtering
int32_t rral_remove_by_prominence(uint16_t* peaks, uint16_t peaksLen, uint16_t peakCount, int32_t* data, uint16_t dataLen, int32_t pkProm, int32_t pkHeval, int32_t pkWidth) {
	if (peakCount == 0) { return RRAL_ERR; }
	if (dataLen == 0) { return RRAL_ERR; }
	if (peakCount >= peaksLen) { peakCount = peaksLen - 1; }

	// The last bit of a peaks index is used as a flag for whether the peak is to be removed after this routine is complete
	// This limits the maximum input signal length that this algorithm can work on to a maximum of 32767 samples.

	// Iterate through all peaks
	for (int32_t i = 0; i < peakCount; i++) {
		uint16_t peak = peaks[i] &= 0x7FFF;
		uint16_t neighbour = 0;

		// Find the left prominence
		int32_t leftProminence = -1;
		for (int32_t j = i - 1; j >= 0; j--) {
			neighbour = peaks[j] & 0x7FFF;
			if (neighbour == 0) { continue; }
			if (data[neighbour] > data[peak]) {
				leftProminence = data[peak] - rral_arr_min(data, dataLen, neighbour, peak);
				break;
			}
		}

		// If there weren't any peaks taller than this one to the left of it then the prominence value is the lowest valley to this peak's left
		if (leftProminence == -1) {
			leftProminence = data[peak] - rral_arr_min(data, dataLen, 0, peak);
		}

		// If prominence is smaller than the minimum, remove this peak
		if (leftProminence < pkProm) {
			peaks[i] |= 0x8000;
			continue;
		}

		// Find the right prominence
		int32_t rightProminence = -1;
		for (int32_t j = i + 1; j < peakCount; j++) {
			neighbour = peaks[j] & 0x7FFF;
			if (neighbour == 0) { continue; }
			if (data[neighbour] > data[peak]) {
				rightProminence = data[peak] - rral_arr_min(data, dataLen, peak, neighbour);
				break;
			}
		}

		// If there weren't any peaks taller than this one to the right of it then the prominence value is the lowest valley to this peak's right
		if (rightProminence == -1) {
			rightProminence = data[peak] - rral_arr_min(data, dataLen, peak, dataLen - 1);
		}

		// If prominence is smaller than the minimum, remove this peak
		if (rightProminence < pkProm) {
			peaks[i] |= 0x8000;
			continue;
		}

		// The smallest value is this peak's prominence; keep it in rightProminence to save 4 bytes.
		if (leftProminence < rightProminence) {
			rightProminence = leftProminence;
		}

		// If prominence is larger than the maximum, remove this peak
		// if(rightProminence > pkPromMax){
		// 	peaks[i] = peaks[i] |= 0x8000;
		// 	continue;
		// }

		// Check this peak's width
		int32_t evaluationHeight = data[peak] - rral_mul32(rightProminence, pkHeval);
		int32_t width = 0;

		// Width on the left
		for (int32_t j = peak; j >= 0; j--) {
			if (j >= dataLen) { break; }
			if (data[j] <= evaluationHeight) {
				width += (peak - j);
				break;
			}
		}

		// Width on the right
		for (int32_t j = peak; j < dataLen; j++) {
			if (j < 0) { break; }
			if (data[j] <= evaluationHeight) {
				width += (j - peak);
				break;
			}
		}

		// Mark peak for removal is width is smaller than pkWidth
		if (rral_int_to_fixed(width) < pkWidth) {
			peaks[i] |= 0x8000;
			continue;
		}
	}

	// Remove marked peaks
	for (uint16_t i = 0; i < peakCount; i++) {
		if ((peaks[i] & 0x8000) != 0) { peaks[i] = 0; }
	}

	return RRAL_OK;
}

// Filter peaks based on proximity // TODO: mark as static?
int32_t rral_remove_close_peaks(uint16_t* peaks, uint16_t peaksLen, int16_t peakCount, int32_t* data, uint16_t dataLen, int32_t dist){
	if(peakCount == 0){ return RRAL_ERR;}
	if(peakCount >= peaksLen){peakCount = peaksLen-1;}

	for(uint16_t i = 0; i < peakCount; i++){
		if(peaks[i] == 0){ continue; }
		for(uint16_t j = i+1; j < peakCount; j++){
			if(peaks[j] == 0){ continue; }
			if(peaks[j]-peaks[i] >= dist){
				break;
			}else{
				if(data[peaks[i]] >= data[peaks[j]]){
					peaks[j] = 0;
				}else{
					peaks[i] = 0;
					break;
				}
			}
		}
	}
	return RRAL_OK;
}

// find_peaks algorithm
int32_t rral_find_peaks(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, Params_findpeaks_t* params) {
	if(data == NULL){return RRAL_ERR;}
	if(peaks == NULL){return RRAL_ERR;}
	if(params == NULL){return RRAL_ERR;}

	// Initialize peak array
	for (uint16_t i = 0; i < peaksLen; i++) {
		peaks[i] = 0;
	}

	uint8_t ret = 0;

	// Find all local maxima
	// uint16_t peakCount = rral_local_maxima(peaks, peaksLen, data, dataLen, pkOrder);
	uint16_t peakCount = rral_local_maxima_fast(peaks, peaksLen, data, dataLen);

	// Filter peaks based on their prominence and width
	ret = rral_remove_by_prominence(peaks, peaksLen, peakCount, data, dataLen, params->pkProm, params->pkHeval, params->pkWidth);
	if (ret == RRAL_ERR) { return RRAL_ERR; }

	// Remove peaks that are too close together(lowest ones removed first)
	ret = rral_remove_close_peaks(peaks, peaksLen, peakCount, data, dataLen, params->pkProxim);
	if (ret == RRAL_ERR) { return RRAL_ERR; }

	// Shift all peaks left to remove any blank space in the buffer
	peakCount = 0;
	uint16_t empty = 0;
	for (uint16_t i = 0; i < peaksLen; i++) {
		if (peaks[i] == 0) { empty++; }
		else if(i - empty != i){ peaks[i - empty] = peaks[i]; peaks[i] = 0; peakCount++; }
		else{ peakCount++; }
	}

	return peakCount;
}

// TERMA algorithm
int32_t rral_terma(int32_t* data, uint16_t dataLen, int32_t bCoef){
	if(data == NULL){return RRAL_ERR;}
	if(EVENT_BUF_SIZE >= dataLen){return RRAL_ERR;}
	if(CYCLE_BUF_SIZE >= dataLen){return RRAL_ERR;}
	if(dataLen <= 0){return RRAL_ERR;}
	
	// Event circular buffer
	int32_t circBufEvent[EVENT_BUF_SIZE];
	uint16_t eventWptr = 0;
	uint16_t eventBufMid = EVENT_BUF_SIZE/2;

	// Cycle circular buffer
	int32_t circBufCycle[CYCLE_BUF_SIZE];
	uint16_t cycleWptr = 0;
	uint16_t cycleBufMid = CYCLE_BUF_SIZE/2;

	// Variables for efficient buffer sum recalculation
	int32_t evOldest = circBufEvent[0];
	int32_t evNewest = circBufEvent[EVENT_BUF_SIZE-1];
	int64_t eventSum = 0;
	int32_t cyOldest = circBufCycle[0];
	int32_t cyNewest = circBufCycle[CYCLE_BUF_SIZE-1];
	int64_t cycleSum = 0;

	// Starting state of the event buffer
	for(uint16_t i = 0; i < EVENT_BUF_SIZE; i++){
		if(i < eventBufMid){
			circBufEvent[i] = 0;
		}else{
			circBufEvent[i] = data[i-eventBufMid];
			eventSum += circBufEvent[i];
		}
	}

	// Starting state of the cycle buffer
	for(uint16_t i = 0; i < CYCLE_BUF_SIZE; i++){
		if(i < cycleBufMid){
			circBufCycle[i] = 0;
		}else{
			circBufCycle[i] = data[i-cycleBufMid];
			cycleSum += circBufCycle[i];
		}
	}

	// Data mean
	int64_t mean = 0;
	for(uint16_t i = 0; i < dataLen; i++){
		mean += data[i];
	}
	mean = mean / dataLen;

	int32_t aCoef = rral_mul32(bCoef, (int32_t)mean);

	// Main iteration
	for(uint16_t i = 0; i < dataLen; i++){
		data[i] = (int64_t)(eventSum / EVENT_BUF_SIZE) - ((cycleSum / CYCLE_BUF_SIZE) + aCoef); // terma_corrected
		// data[i] = (int64_t)(eventSum / EVENT_BUF_SIZE) * 2 + ((cycleSum / CYCLE_BUF_SIZE) + aCoef); // terma_star

		// Update event buffer
		if(i+eventBufMid >= dataLen){
			// Pad with 0s if not enough data left
			evOldest = circBufEvent[eventWptr];
			circBufEvent[eventWptr] = 0;
			evNewest = 0;
		}else{
			evOldest = circBufEvent[eventWptr];
			circBufEvent[eventWptr] = data[i+eventBufMid];
			evNewest = circBufEvent[eventWptr];
		}

		// Update cycle buffer
		if(i+cycleBufMid >= dataLen){
			// Pad with 0s if not enough data left
			cyOldest = circBufCycle[cycleWptr];
			circBufCycle[cycleWptr] = 0;
			cyNewest = 0;
		}else{
			cyOldest = circBufCycle[cycleWptr];
			circBufCycle[cycleWptr] = data[i+cycleBufMid];
			cyNewest = circBufCycle[cycleWptr];
		}

		// Wrap condition
		eventWptr = (eventWptr + 1) % EVENT_BUF_SIZE;
		cycleWptr = (cycleWptr + 1) % CYCLE_BUF_SIZE;

		// Recalculate circular buffer sums
		eventSum = eventSum - evOldest + evNewest;
		cycleSum = cycleSum - cyOldest + cyNewest;
	}

	return RRAL_OK;
}

// SRMAC algorithm
int32_t rral_srmac(int32_t* data, uint16_t dataLen, Params_srmac_t params, State_srmac_t *state) {
	if(data == NULL){return RRAL_ERR;}
	if(state == NULL){return RRAL_ERR;}
	if (dataLen <= 0) { return RRAL_ERR; }

	// State variables
	int64_t fast = 0;
	int64_t slow = 0;
	int64_t cross = 0;

	if (state != NULL) { 
		// Dereference state if available
		fast = state->prevFast;
		slow = state->prevSlow;
		cross = state->prevCross;
	}

	// Alt coefficients
	int32_t coefFast_b = rral_int_to_fixed(1) - params.coefFast;
	int32_t coefSlow_b = rral_int_to_fixed(1) - params.coefSlow;
	int32_t coefCros_b = rral_int_to_fixed(1) - params.coefCross;

	// Intermediary data
	int32_t curr = 0;

	// Calculations
	for (uint16_t i = 0; i < dataLen; i++) {
		curr = data[i];
		fast = (((int64_t)curr * params.coefFast) >> FIXED_POINT) + (((int64_t)fast * coefFast_b) >> FIXED_POINT);
		slow = (((int64_t)curr * params.coefSlow) >> FIXED_POINT) + (((int64_t)slow * coefSlow_b) >> FIXED_POINT);
		cross = ((((int64_t)fast - slow) * params.coefCross) >> FIXED_POINT) + (((int64_t)cross * coefCros_b) >> FIXED_POINT);
		if(cross > 0x7FFFFFFF){cross = 0x7FFFFFFF;} 	// clamp ovf
		// if(cross < -0x80000000){cross = -0x80000000;} 	// ...
		data[i] = cross;
	}

	if (state != NULL) {
		// Save state
		state->prevFast = fast;
		state->prevSlow = slow;
		state->prevCross = cross;
	}

	return RRAL_OK;
}

// Find peaks in zero-crossing data
int32_t rral_zero_crossing(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th){
	if (dataLen <= 0) { return RRAL_ERR;}
	if (data == NULL) { return RRAL_ERR;}
	if (peaks == NULL) { return RRAL_ERR;}

	uint16_t count = 0;
	int32_t positive = 0;
	int32_t distToPeak = 0;
	int32_t maximum = 0;

	for(uint16_t i = 0; i < dataLen; i++){
		if(data[i] > th){ // Positive edge
			positive++;
			distToPeak++;
			if(maximum < data[i]){
				maximum = data[i];
				distToPeak = 0;
			}
		}else{ // Negative edge
			if(positive >= width){ // Real peak detected
				peaks[count] = i-distToPeak-1;
				count++;
				if(count >= peaksLen){count = 0;}
			}
			positive = 0;
			distToPeak = 0;
			maximum = 0;
		}
	}
	return count;
}

// Find peaks in zero-crossing data with peak maxima taken from raw unfiltered data
int32_t rral_zero_crossing_raw(int32_t* data, int32_t* rawData, uint16_t dataLen, uint16_t* peaks, uint16_t peaksLen, int32_t width, int32_t th){
	if (dataLen <= 0) { return RRAL_ERR;}
	if (data == NULL) { return RRAL_ERR;}
	if (peaks == NULL) { return RRAL_ERR;}
	if (rawData == NULL) { return RRAL_ERR;}

	uint16_t count = 0;
	int32_t positive = 0;
	int32_t distToPeak = 0;
	int32_t maximum = 0;

	for(uint16_t i = 0; i < dataLen; i++){
		if(data[i] > th){ // Positive edge
			positive++;
			distToPeak++;
			if(maximum < rawData[i]){
				maximum = rawData[i];
				distToPeak = 0;
			}
		}else{ // Negative edge
			if(positive >= width){ // Real peak detected
				peaks[count] = i-distToPeak-1;
				count++;
				if(count >= peaksLen){count = 0;}
			}
			positive = 0;
			distToPeak = 0;
			maximum = 0;
		}
	}
	return count;
}

// Neighbour running covariance SQI method - returns the window Pearson coefficient
int32_t rral_get_sqi_lite(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount, int32_t mindist, int32_t maxdist){
	if(data == NULL){return RRAL_ERR;}
	if(peaks == NULL){return RRAL_ERR;}
	if(peakCount < 2){ return 0; }

	// Calculate the average distance between peaks
	int32_t dist = (peaks[peakCount-1]-peaks[0]) / (peakCount - 1);
	if(dist % 2 != 0){ dist++; }

	int64_t coefs = 0;
	int32_t coefsN = 0;
	int32_t inv = 0;

	// For each pair of peaks calculate their running covariance using the online covariance algorithm as described in https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online
	for(int32_t i = 0; i < peakCount-1; i++){

		int32_t halfdist = dist / 2;

		if(peaks[i] < halfdist){ inv++; continue;}
		if(dataLen-peaks[i+1] < halfdist){ inv++; continue;}
		if(halfdist < SAMPLE_RATE / 4){ inv++; continue; }

		int32_t* fst = &data[peaks[i]];
		int32_t* snd = &data[peaks[i+1]];
		int64_t dx = 0;
		int64_t dx2 = 0;
		int64_t dy = 0;
		int64_t dy2 = 0;
		int64_t meanx = 0;
		int64_t meany = 0;
		int64_t cov = 0;
		int32_t n = 0;
		int64_t varx = 0;
		int64_t vary = 0;

		// First pass to get the mean values
		for(int32_t j = -halfdist; j < halfdist; j++){
			n++;
			meanx += fst[j];
			meany += snd[j];
		}
		meanx /= n;
		meany /= n;
		n = 0;

		// Second pass for live calculation
		int8_t scale = 0;
		for(int32_t j = -halfdist; j < halfdist; j++){
			n++;

			dx = fst[j] - meanx;
			dy = snd[j] - meany;

			dx >>= scale;
			dy >>= scale;

			// cov
			cov += (dx * dy) >> scale;

			// var
			varx += (dx * dx) >> scale;
			vary += (dy * dy) >> scale;
		}

		// // Integer version
		cov /= (n-1);
		varx /= (n-1);
		vary /= (n-1);

		if(varx != 0 && vary != 0){
			// cov = (cov << FIXED_POINT) / ((varx+vary)/2); // if we assume that varx ~= vary then sqrt(varx*vary) ~= (varx+vary)/2 (very rough approx. that saves time)
			// coefs += cov;
			// coefsN++;

			// Accurate coefficient using double math
			int64_t prod = 0;
			double d_prod = (double)varx / (float)(1 << FIXED_POINT);
			d_prod *= (double)vary / (float)(1 << FIXED_POINT);
			d_prod = sqrt(d_prod);
			prod = (int64_t)(d_prod * (double)(1 << FIXED_POINT) + (d_prod >= 0 ? 0.5 : -0.5));
			cov = (cov << FIXED_POINT) / prod;
			coefs += cov;
			coefsN++;
		}
	}

	if(inv >= peakCount-1){
		return 0;
	}

	if(coefsN == 0){
		return RRAL_ERR;
	}

	return coefs / coefsN;
}

// Full SQI method - returns the window Pearson coefficient
int32_t rral_get_sqi_full(int32_t* data, uint16_t dataLen, uint16_t* peaks, uint16_t peakCount){
	if(data == NULL){return RRAL_ERR;}
	if(peaks == NULL){return RRAL_ERR;}
	
	// Calculate the median distance between peaks
	int32_t dist = 0;
	if(peakCount < 2){ return 0; }
	else{
		dist = (peaks[peakCount-1]-peaks[0]) / (peakCount - 1);
	}
	if(dist % 2 != 0){ dist++; }

	// Are < 15% of the breath durations > 1.5 or < 0.5 times the median breath duration?
	int32_t pn = 0;
	for(uint16_t i = 0; i < peakCount-1; i++){
		int32_t tdist = peaks[i+1] - peaks[i];
		if(tdist > ((dist * 3) / 2) || tdist < (dist / 2)){
			pn++;
		}
	}
	// e.g if pn / peakcount > 0.15 return 0

	// Is the normalized standard deviation of breath durations < 0.25



	// Calculate the correlation coefficient

	// Aux
	int32_t pattern[SAMPLE_RATE * DATA_WINDOW/2];
	int32_t halfdist = dist/2;
	int16_t fst = -1;
	pn = 0;

	// Zero out the pattern buffer
	for(uint16_t i = 0; i < dist; i++){
		pattern[i] = 0;
	}

	// Find first usable peak // Used when not using zero-padding
	for(uint16_t i = 0; i < peakCount; i++){
		if(peaks[i] < halfdist || halfdist > dataLen-peaks[i]){ continue; }
		if(fst == -1){fst = i;}
		pn++;
	}

	// These are the initial variables for the zero-padding version
	fst = 0;
	pn = peakCount;
	uint16_t idx = 0;

	// Calculate the average pattern
	for(uint16_t i = 0; i < dist; i++){
		int64_t temp = 0;
		for(uint16_t j = fst; j < pn; j++){
			// temp += data[peaks[fst+j]-halfdist+i]; // Previous version - I think this is wrong but not 100% sure so leaving it in

			// Zero-padding when peaks are too close to the edge
			idx = peaks[j]-halfdist+i;
			if(idx >= 0 && idx < dataLen){
				temp += data[peaks[j]-halfdist+i];
			}
		}
		pattern[i] = temp / pn;
	}

	// For each peak calculate the Pearson correlation coefficient with the average pattern
	// OPTIMIZATION: the mean of the pattern does not need to be recalculated every time; neither does vary; get rid of double math;
	// OPTIMIZATION: zero-padding increases computational load but can be mitigated by performing a check before each for loop and applying a loop without an if statement if peak is within the bounds
	int64_t coefs = 0;
	int32_t validCoefs = 0;
	int32_t inv = 0;
	for(int32_t i = fst; i < pn; i++){

		int32_t* fst = &data[peaks[i]];
		int32_t* snd = &pattern[halfdist];
		int64_t dx = 0;
		int64_t dy = 0;
		int64_t meanx = 0;
		int64_t meany = 0;
		int64_t cov = 0;
		int64_t varx = 0;
		int64_t vary = 0;

		// First pass to get the mean values
		for(int32_t j = -halfdist; j < halfdist; j++){
			// Zero-padding when peaks are too close to the edge
			if(peaks[i]+j >= 0 && peaks[i]+j < dataLen){
				meanx += fst[j];
			}
			meany += snd[j];
		}
		meanx /= dist;
		meany /= dist;

		// Second pass for live calculation
		int8_t scale = 0;
		for(int32_t j = -halfdist; j < halfdist; j++){
			// Zero-padding when peaks are too close to the edge
			if(peaks[i]+j >= 0 && peaks[i]+j < dataLen){
				dx = fst[j] - meanx;
			}
			dy = snd[j] - meany;

			dx >>= scale;
			dy >>= scale;

			// cov
			cov += (dx * dy) >> scale;

			// var
			varx += (dx * dx) >> scale;
			vary += (dy * dy) >> scale;
		}

		// // Integer version
		cov /= (dist-1);
		varx /= (dist-1);
		vary /= (dist-1);

		if(varx != 0 || vary != 0){
			// Very rough coefficient approximation using just integers
			// cov = (cov << FIXED_POINT) / ((varx+vary)/2); // if we assume that varx ~= vary then sqrt(varx*vary) ~= (varx+vary)/2

			// Accurate coefficient using double math
			int64_t prod = 0;
			double d_prod = (double)varx / (float)(1 << FIXED_POINT);
			d_prod *= (double)vary / (float)(1 << FIXED_POINT);
			d_prod = sqrt(d_prod);
			prod = (int64_t)(d_prod * (double)(1 << FIXED_POINT) + (d_prod >= 0 ? 0.5 : -0.5));
			cov = (cov << FIXED_POINT) / prod;
			coefs += cov;
			validCoefs++;
		}
	}

	if(inv >= peakCount-1){
		return 0;
	}

	if(validCoefs == 0){
		return RRAL_ERR;
	}

	return coefs / validCoefs;
}

// Zero-init single parameter state structs
int32_t rral_single_param_state_init(State_joint_t* stateJ, State_peaks_t* state){
	if(stateJ == NULL){return RRAL_ERR;}
	if(state == NULL){return RRAL_ERR;}

	// Init data buffers
	for(int32_t i = 0; i < MAIN_BUF_SIZE; i++){
		state->data[i] = 0;
	}

	// Init peak buffer
	for(int32_t i = 0; i < MAX_PEAK_COUNT; i++){
		state->peaks[i] = 0;
	}

	// Init SRMAC states
	state->srmacState.prevFast = 0;
	state->srmacState.prevSlow = 0;
	state->srmacState.prevCross = 0;

	// Init filter state buffers
	for(int32_t i = 0; i < 3*SOS_SECT; i++){
		state->hpfx[i] = 0;
		state->hpfy[i] = 0;
		state->lpfx[i] = 0;
		state->lpfy[i] = 0;
	}

	// Init running values
	state->deltaSampleCount = 0;
	stateJ->deltaSampleCount = 0;
	stateJ->totalSampleCount = 0;
	stateJ->runningMeanSum = 0;
	stateJ->runningMean = 0;

	return RRAL_OK;
}

// Zero-init joint hr/rr state structs
int32_t rral_joint_hr_rr_state_init(State_joint_t* stateJ, State_peaks_t* stateHR, State_peaks_t* stateRR){
	if(stateJ == NULL){return RRAL_ERR;}
	if(stateHR == NULL){return RRAL_ERR;}
	if(stateRR == NULL){return RRAL_ERR;}

	// Init data buffers
	for(int32_t i = 0; i < MAIN_BUF_SIZE; i++){
		stateHR->data[i] = 0;
		stateRR->data[i] = 0;
	}

	// Init peak buffer
	for(int32_t i = 0; i < MAX_PEAK_COUNT; i++){
		stateHR->peaks[i] = 0;
		stateRR->peaks[i] = 0;
	}

	// Init SRMAC states
	stateHR->srmacState.prevFast = 0;
	stateHR->srmacState.prevSlow = 0;
	stateHR->srmacState.prevCross = 0;
	stateRR->srmacState.prevFast = 0;
	stateRR->srmacState.prevSlow = 0;
	stateRR->srmacState.prevCross = 0;

	// Init filter state buffers
	for(int32_t i = 0; i < 3*SOS_SECT; i++){
		stateHR->hpfx[i] = 0;
		stateHR->hpfy[i] = 0;
		stateHR->lpfx[i] = 0;
		stateHR->lpfy[i] = 0;
		stateRR->hpfx[i] = 0;
		stateRR->hpfy[i] = 0;
		stateRR->lpfx[i] = 0;
		stateRR->lpfy[i] = 0;
	}

	// Init running values
	stateJ->deltaSampleCount = 0;
	stateHR->deltaSampleCount = 0;
	stateRR->deltaSampleCount = 0;
	stateJ->totalSampleCount = 0;
	stateJ->runningMeanSum = 0;
	stateJ->runningMean = 0;

	return RRAL_OK;
}

// Realtime peak detection block - buffers need to be managed externally - this function only filters the last RRAL_DELTA_SAMPLES samples and detects peaks
int32_t rral_realtime_peaks(Params_peaks_t* params, State_peaks_t* state){
	if(params == NULL){return RRAL_ERR;}
	if(state == NULL){return RRAL_ERR;}

	// Clear peak buffer
	for(uint16_t i = 0; i < MAX_PEAK_COUNT; i++){ 
		state->peaks[i] = 0;
	}

	// HR bandpass filter
	rral_sosfilt(&state->data[MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES], RRAL_DELTA_SAMPLES, params->hpfsos, SOS_SECT, state->hpfx, state->hpfy);
	rral_sosfilt(&state->data[MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES], RRAL_DELTA_SAMPLES, params->lpfsos, SOS_SECT, state->lpfx, state->lpfy);

	// Find heartbeats using SRMAC
	rral_srmac(&state->data[MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES], RRAL_DELTA_SAMPLES, params->srmacParams, &state->srmacState);

	// Inline zero-crossing detection
	state->peaksN = rral_zero_crossing(&state->data[0], MAIN_BUF_SIZE, state->peaks, MAX_PEAK_COUNT, params->srmacWidth, 0);

	return RRAL_OK;
}

// Shift the algorithm data buffer
int32_t rral_realtime_peaks_shift(State_peaks_t* state){
	if(state == NULL){return RRAL_ERR;}

	// Shift all samples in data back by RRAL_DELTA_SAMPLES
	for(uint16_t i = RRAL_DELTA_SAMPLES; i < MAIN_BUF_SIZE; i++){
		state->data[i - RRAL_DELTA_SAMPLES] = state->data[i];
	}

	return RRAL_OK;
}

// Single parameter estimation wrapper
int32_t rral_single_param(int16_t sample, Params_peaks_t* params, State_peaks_t* state, State_joint_t* stateJ, Results_single_t* out){
	if(params == NULL){return RRAL_ERR;}
	if(state == NULL){return RRAL_ERR;}
	if(stateJ == NULL){return RRAL_ERR;}
	if(out == NULL){return RRAL_ERR;}

	int32_t f_sample = 0;

	// Subtract the running mean from the sample
	if(stateJ->totalSampleCount < RUNNING_MEAN_N){
		stateJ->runningMeanSum += sample;
		stateJ->runningMean = stateJ->runningMeanSum / stateJ->totalSampleCount;
	}else{
		stateJ->runningMeanSum -= stateJ->runningMean;
		stateJ->runningMeanSum += sample;
		stateJ->runningMean = stateJ->runningMeanSum / RUNNING_MEAN_N;
	}

	f_sample = sample - stateJ->runningMean;
	f_sample = rral_int_to_fixed(f_sample);

	// Add the sample to the data buffers
	state->data[(MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES) + stateJ->deltaSampleCount] = f_sample;
	stateJ->deltaSampleCount++;
	stateJ->totalSampleCount++;

	if(stateJ->deltaSampleCount == RRAL_DELTA_SAMPLES){
		stateJ->deltaSampleCount = 0;

		Results_single_t results = {0, 0, 0};

		results.time = (stateJ->totalSampleCount - RRAL_DELTA_SAMPLES) / SAMPLE_RATE;

		// -- Rate calculation -- //

		uint16_t count = 0;
		if(rral_realtime_peaks(params, state) == RRAL_OK){
			count = state->peaksN;
		}

		// Calculate the rate over the entire buffer window
		if(count > 2){ // If possible skip the left-most peak as it will often have been imprecisely detected
			results.rate = rral_int_to_fixed(state->peaks[count - 1] - state->peaks[1]) / (count - 2);
			results.rate = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.rate);
			results.sqi = rral_get_sqi_full(&state->data[0], MAIN_BUF_SIZE, state->peaks, count);
		}else if(count > 1){
			results.rate = rral_int_to_fixed(state->peaks[count - 1] - state->peaks[0]) / (count - 1);
			results.rate = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.rate);
			results.sqi = rral_get_sqi_full(&state->data[0], MAIN_BUF_SIZE, state->peaks, count);
		}else{
			results.rate = 0;
			results.sqi = 0;
		}

		// Perform post-estimate shift
		rral_realtime_peaks_shift(state);

		// Return results
		*out = results;

		return RRAL_RDY;
	}

	return RRAL_OK;
}

// Realtime joint HR and RR algorithm
int32_t rral_joint_hr_rr(int16_t sample, Params_peaks_t* paramsHR, Params_peaks_t* paramsRR, State_peaks_t* stateHR, State_peaks_t* stateRR, State_joint_t* stateJ, Results_joint_t* out){
	if(paramsHR == NULL){return RRAL_ERR;}
	if(paramsRR == NULL){return RRAL_ERR;}
	if(stateHR == NULL){return RRAL_ERR;}
	if(stateRR == NULL){return RRAL_ERR;}
	if(stateJ == NULL){return RRAL_ERR;}
	if(out == NULL){return RRAL_ERR;}

	int32_t f_sample = 0;

	// Subtract the running mean from the sample
	if(stateJ->totalSampleCount < RUNNING_MEAN_N){
		stateJ->runningMeanSum += sample;
		stateJ->runningMean = stateJ->runningMeanSum / stateJ->totalSampleCount;
	}else{
		stateJ->runningMeanSum -= stateJ->runningMean;
		stateJ->runningMeanSum += sample;
		stateJ->runningMean = stateJ->runningMeanSum / RUNNING_MEAN_N;
	}

	f_sample = sample - stateJ->runningMean;
	f_sample = rral_int_to_fixed(f_sample);

	// Add the sample to the data buffers
	stateHR->data[(MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES) + stateJ->deltaSampleCount] = f_sample;
	stateRR->data[(MAIN_BUF_SIZE - RRAL_DELTA_SAMPLES) + stateJ->deltaSampleCount] = f_sample;
	stateJ->deltaSampleCount++;
	stateJ->totalSampleCount++;

	if(stateJ->deltaSampleCount == RRAL_DELTA_SAMPLES){
		stateJ->deltaSampleCount = 0;

		Results_joint_t results = {0, 0, 0, 0, 0};

		results.time = (stateJ->totalSampleCount - RRAL_DELTA_SAMPLES) / SAMPLE_RATE;

		// -- HR routine -- //

		uint16_t beatcount = 0;
		if(rral_realtime_peaks(paramsHR, stateHR) == RRAL_OK){
			beatcount = stateHR->peaksN;
		}

		// Calculate the HR over the entire buffer window
		if(beatcount > 2){ // If possible skip the left-most heartbeat as it will often have been imprecisely detected due to edge effects
			results.hr = rral_int_to_fixed(stateHR->peaks[beatcount - 1] - stateHR->peaks[1]) / (beatcount - 2);
			results.hr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.hr);
			// results.hr = rral_int_to_fixed(stateHR->peaks[beatcount - 1] - stateHR->peaks[beatcount - 2]);
			// results.hr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.hr);
			results.sqihr = rral_get_sqi_full(&stateHR->data[0], MAIN_BUF_SIZE, stateHR->peaks, beatcount);
		}else if(beatcount > 1){
			results.hr = rral_int_to_fixed(stateHR->peaks[beatcount - 1] - stateHR->peaks[0]) / (beatcount - 1);
			results.hr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.hr);
			// results.hr = rral_int_to_fixed(stateHR->peaks[beatcount - 1] - stateHR->peaks[beatcount - 2]);
			// results.hr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.hr);
			results.sqihr = rral_get_sqi_full(&stateHR->data[0], MAIN_BUF_SIZE, stateHR->peaks, beatcount);
		}else{
			results.hr = 0;
			results.sqihr = 0;
		}

		// Perform post-estimate shift
		rral_realtime_peaks_shift(stateHR);

		// --- HR end ----- //

		// -- RR routine -- //

		uint16_t breathcount = 0;
		if(rral_realtime_peaks(paramsRR, stateRR) == RRAL_OK){
			breathcount = stateRR->peaksN;
		}

		// Calculate the RR over the entire buffer window
		if(breathcount > 2){ // If possible skip the left-most breath as it will often have been imprecisely detected due to edge effects
			results.rr = rral_int_to_fixed(stateRR->peaks[breathcount - 1] - stateRR->peaks[1]) / (breathcount - 2);
			results.rr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.rr);
			results.sqirr = rral_get_sqi_full(&stateRR->data[0], MAIN_BUF_SIZE, stateRR->peaks, breathcount);
		}else if(breathcount > 1){
			results.rr = rral_int_to_fixed(stateRR->peaks[breathcount - 1] - stateRR->peaks[0]) / (breathcount - 1);
			results.rr = rral_div32(rral_int_to_fixed(60 * SAMPLE_RATE), results.rr);
			results.sqirr = rral_get_sqi_full(&stateRR->data[0], MAIN_BUF_SIZE, stateRR->peaks, breathcount);
		}else{
			results.rr = 0;
			results.sqirr = 0;
		}

		// Perform post-estimate shift
		rral_realtime_peaks_shift(stateRR);

		// --- RR end ----- //

		// Return results
		*out = results;

		return RRAL_RDY;
	}

	return RRAL_OK;
}

#if RRAL_CWT_ENABLED == 1
// Mexh wavelet
float rral_mexh(float t){
	// 0.867 is the Mexh wavelet normalization constant // ref: https://ieeexplore.ieee.org/document/9941083
	return 0.867*powf(2.72, -(t*t/2))*(1-(t*t));
}

// Gaus2 wavelet
float rral_gaus2(float t){
	// 1.031 is the Gaus2 wavelet normalization constant // ref: https://ieeexplore.ieee.org/document/9941083
	return 1.031*powf(2.72, -(t*t))*(1-2*(t*t));
}

// Scaled version of the Mexh wavelet
float rral_mexh_scaled(float t, float scale){
	if(scale == 0){return 0;}
	return (1/sqrtf(scale))*rral_mexh(t/scale);
}

// Scaled version of the Gaus2 wavelet
float rral_gaus2_scaled(float t, float scale){
	if(scale == 0){return 0;}
	return (1/sqrtf(scale))*rral_gaus2(t/scale);
}

// Scaled wavelet wrapper function
float rral_wavelet_scaled(float t, float scale){
	if(scale == 0){return 0;}
	#if RRAL_CWT_WAVELET == 1
	return (1/sqrtf(scale))*rral_mexh(t/scale);
	#else
	return (1/sqrtf(scale))*rral_gaus2(t/scale);
	#endif
}

// CWT implementation using almost entirely floating point calculations to take advantage of the existing dsps_conv.h signal processing methods
int32_t rral_cwt(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table){
	if(data == NULL){return RRAL_ERR;}
	if(fft_table == NULL){return RRAL_ERR;}

	static float scales[CWT_SCALES_SIZE];

	static float signal[CWT_BUF_SIZE];
	static float kernel[CWT_BUF_SIZE];
	static float result[CWT_BUF_SIZE];

	// Calculate scales
	float scales_step = (maxFreq-minFreq)/CWT_SCALES_SIZE;
	for(uint16_t i = 0; i < CWT_SCALES_SIZE; i++){
		scales[i] = 0.25/((minFreq+i*scales_step)/SAMPLE_RATE);
	}

	// Instantiate signal array
	for(uint16_t i = 0; i < CWT_BUF_SIZE; i++){
		signal[i] = 0;
	}

	// Find maximum amount of signal data that the signal buffer can contain(everything after has to be zero-padded to avoid circular convolution)
	uint16_t cwtDataLen = CWT_BUF_SIZE - KERNEL_SIZE;

	// Move data to signal array
	for(uint16_t i = 0; i < dataLen; i++){
		if(i*2+1 < cwtDataLen){
			signal[i*2] = rral_fixed_to_float(data[i]);
			data[i] = 0; // Zero out the initial data array
		}else{
			data[i] = 0; // Zero out the initial data array
		}
	}	

	// Precompute the FFT of the input signal
	dsps_fft2r_fc32_ansi_(signal, CWT_BUF_SIZE/2, fft_table);
	dsps_bit_rev_fc32_ansi(signal, CWT_BUF_SIZE/2);

	// Iterate through each scale
	for(uint8_t k = 0; k < CWT_SCALES_SIZE; k++){

		// Select scale
		float scale = scales[k];

		// Zero out the kernel buffer
		for(uint16_t i = 0; i < CWT_BUF_SIZE; i++){
			kernel[i] = 0;
		}

		// Fill all even indices of the kernel buffer with real values of the chosen wavelet
		for(int16_t i = 0; i < KERNEL_SIZE; i++){
			if(i*2 < CWT_BUF_SIZE){
				kernel[i*2] = rral_wavelet_scaled(((float)i)-(KERNEL_SIZE/2), scale);
			}
		}

		// Compute the FFT of the kernel
		dsps_fft2r_fc32_ansi_(kernel, CWT_BUF_SIZE/2, fft_table);
		dsps_bit_rev_fc32_ansi(kernel, CWT_BUF_SIZE/2);

		// Perform convolution and return the complex conjugate so the forward DFT can be used as an IDFT 
		// (NOTE: swapping real and imaginary values once before the FFT call and once after lets you skip the conjugation step and achieves the same results however this is slightly slower)
		for(uint16_t i = 0; (i*2+1) < CWT_BUF_SIZE; i++){
			// Convolution
			float rv = signal[i*2]*kernel[i*2] - signal[i*2+1]*kernel[i*2+1];
			float iv = signal[i*2]*kernel[i*2+1] + signal[i*2+1]*kernel[i*2];

			// Complex conjugate of the result
			result[i*2] = rv;
			result[i*2+1] = -iv;
		}

		// Compute the reverse FFT
		dsps_fft2r_fc32_ansi_(result, CWT_BUF_SIZE/2, fft_table);
		dsps_bit_rev_fc32_ansi(result, CWT_BUF_SIZE/2);

		// Add real results to output
		int32_t N = CWT_BUF_SIZE/2;
		for(uint16_t i = 0; i < dataLen; i++){
			if((KERNEL_SIZE+i*2) < CWT_BUF_SIZE){ // Account for the fact that result data is shifted forward by KERNEL_SIZE/2 (and -2 ?); since the array contains real/imag pairs the actual index shift is simply KERNEL_SIZE
				data[i] += rral_float_to_fixed((result[KERNEL_SIZE+i*2])/N);
			}else{
				data[i] = 0;
			}
		}
	}

	// Calculate the average value for each point in time
	int32_t d = rral_int_to_fixed(CWT_SCALES_SIZE);
	for(uint16_t i = 0; i < dataLen; i++){
		if(data[i] != 0){
			data[i] = rral_div32(data[i], d);
		}
	}

	return RRAL_OK;
}

// Overlap-add CWT to decrease the memory constraints of the algorithm by an order of magnitude
// Could be optimized for speed by storing kernels in memory which would increase memory consumption but increase execution speed massively
int32_t rral_cwt_oa(int32_t* data, uint16_t dataLen, float minFreq, float maxFreq, float* fft_table) {
	if(data == NULL){return RRAL_ERR;}
	if(fft_table == NULL){return RRAL_ERR;}

	float scales[CWT_SCALES_SIZE];

	float window[CWT_OA_SIZE];
	float kernel[CWT_OA_SIZE];
	float result[CWT_OA_SIZE];
	int32_t overlap[KERNEL_SIZE];

	int32_t finalOffset = -KERNEL_SIZE / 2; // Account for phase offset caused by convolution

	// Initialize overlap buffer
	for (uint16_t i = 0; i < KERNEL_SIZE; i++) {
		overlap[i] = 0;
	}

	// Calculate scales
	float scales_step = (maxFreq - minFreq) / CWT_SCALES_SIZE;
	for (uint16_t i = 0; i < CWT_SCALES_SIZE; i++) {
		scales[i] = 0.25 / ((minFreq + i * scales_step) / SAMPLE_RATE);
	}

	// Step through data in windows of size KERNEL_SIZE
	for (int32_t step = 0; step < dataLen; step += KERNEL_SIZE) {

		// Zero out window buffer
		for (int32_t j = 0; j < CWT_OA_SIZE; j++) {
			window[j] = 0;
		}

		// Move data to window buffer and in the data buffer zero out whatever has been moved
		for (int32_t j = 0; j < KERNEL_SIZE; j++) {
			if (step + j < dataLen) {
				window[j * 2] = rral_fixed_to_float(data[step + j]);
				data[step + j] = 0;
			}

			// Add the previous overlap region to the zeroed out segment
			if (step + j + finalOffset >= 0) {
				data[step + j + finalOffset] += overlap[j];
			}

			overlap[j] = 0;
		}

		// Precompute the FFT of the signal segment
		dsps_fft2r_fc32_ansi_(window, CWT_OA_SIZE / 2, fft_table);
		dsps_bit_rev_fc32_ansi(window, CWT_OA_SIZE / 2);

		// Iterate through each scale
		for (uint8_t k = 0; k < CWT_SCALES_SIZE; k++) {

			// Select scale
			float scale = scales[k];

			// Zero out the kernel buffer
			for (uint16_t j = 0; j < CWT_OA_SIZE; j++) {
				kernel[j] = 0;
			}

			// Fill all even indices of the kernel buffer with real values of the chosen wavelet
			for (int16_t j = 0; j < KERNEL_SIZE; j++) {
				if (j * 2 < CWT_OA_SIZE) {
					kernel[j * 2] = rral_wavelet_scaled(((float)j) - (KERNEL_SIZE / 2), scale);
				}
			}

			// Compute the FFT of the kernel
			dsps_fft2r_fc32_ansi_(kernel, CWT_OA_SIZE / 2, fft_table);
			dsps_bit_rev_fc32_ansi(kernel, CWT_OA_SIZE / 2);

			// Perform convolution and return the complex conjugate so the forward DFT can be used as an IDFT 
			for (uint16_t j = 0; (j * 2 + 1) < CWT_OA_SIZE; j++) {
				float rv = window[j * 2] * kernel[j * 2] - window[j * 2 + 1] * kernel[j * 2 + 1];
				float iv = window[j * 2] * kernel[j * 2 + 1] + window[j * 2 + 1] * kernel[j * 2];

				// Complex conjugate of the result
				result[j * 2] = rv;
				result[j * 2 + 1] = -iv;
			}

			// Compute the reverse FFT
			dsps_fft2r_fc32_ansi_(result, CWT_OA_SIZE / 2, fft_table);
			dsps_bit_rev_fc32_ansi(result, CWT_OA_SIZE / 2);

			// Add real results to output
			for (uint16_t j = 0; j < KERNEL_SIZE; j++) {
				if (step + j + finalOffset < dataLen && step + j + finalOffset >= 0) {
					data[step + j + finalOffset] += rral_float_to_fixed(result[j * 2] / (CWT_OA_SIZE / 2));
				}
			}

			// Add real overlap to overlap buffer
			for (uint16_t j = 0; j < KERNEL_SIZE; j++) { // shouldn't this be KERNEL_SIZE-1?
				overlap[j] += rral_float_to_fixed(result[(KERNEL_SIZE + j) * 2] / (CWT_OA_SIZE / 2));
			}
		}
	}

	// Calculate the average value for each point in time
	int32_t d = rral_int_to_fixed(CWT_SCALES_SIZE);
	for (uint16_t i = 0; i < dataLen; i++) {
		if (data[i] != 0) {
			data[i] = rral_div32(data[i], d);
		}
	}

	return RRAL_OK;
}

// Kernel precalculation function for the optimized overlap-add CWT algorithm
int32_t rral_cwt_oa_fast_init(float* kernels, int32_t kernelBufLen, float minFreq, float maxFreq, float* fft_table) {
	if (kernels == NULL) { return RRAL_ERR; }
	if (fft_table == NULL) { return RRAL_ERR; }
	if (kernelBufLen < CWT_OA_SIZE * CWT_SCALES_SIZE) { return RRAL_ERR; }

	float scales[CWT_SCALES_SIZE];

	// Calculate scales
	float scales_step = (maxFreq - minFreq) / CWT_SCALES_SIZE;
	for (uint16_t i = 0; i < CWT_SCALES_SIZE; i++) {
		scales[i] = 0.25 / ((minFreq + i * scales_step) / SAMPLE_RATE);
	}

	// Zero out the kernels buffer
	for (uint16_t i = 0; i < CWT_OA_SIZE * CWT_SCALES_SIZE; i++) {
		kernels[i] = 0;
	}

	// Iterate through each scale
	for (uint8_t k = 0; k < CWT_SCALES_SIZE; k++) {

		// Select scale
		float scale = scales[k];

		// Select kernel
		float* kernel = &kernels[CWT_OA_SIZE * k];

		// Zero out the kernel buffer
		for (uint16_t j = 0; j < CWT_OA_SIZE; j++) {
			kernel[j] = 0;
		}

		// Fill all even indices of the kernel buffer with real values of the chosen wavelet
		for (int16_t j = 0; j < KERNEL_SIZE; j++) {
			if (j * 2 < CWT_OA_SIZE) {
				kernel[j * 2] = rral_wavelet_scaled(((float)j) - (KERNEL_SIZE / 2), scale);
			}
		}

		// Compute the FFT of the kernel
		dsps_fft2r_fc32_ansi_(kernel, CWT_OA_SIZE / 2, fft_table);
		dsps_bit_rev_fc32_ansi(kernel, CWT_OA_SIZE / 2);
	}

	return RRAL_OK;
}

// Overlap-add CWT with saved kernels for better performance
int32_t rral_cwt_oa_fast(int32_t* data, uint16_t dataLen, float* fft_table, float* kernels) {
	if(data == NULL){return RRAL_ERR;}
	if(kernels == NULL){return RRAL_ERR;}
	if(fft_table == NULL){return RRAL_ERR;}

	float window[CWT_OA_SIZE];
	float result[CWT_OA_SIZE];
	int32_t overlap[KERNEL_SIZE];

	int32_t finalOffset = -KERNEL_SIZE / 2; // Account for phase offset caused by convolution

	// Initialize overlap buffer
	for (uint16_t i = 0; i < KERNEL_SIZE; i++) {
		overlap[i] = 0;
	}

	// Step through data in windows of size KERNEL_SIZE
	for (int32_t step = 0; step < dataLen; step += KERNEL_SIZE) {

		// Zero out window buffer
		for (int32_t j = 0; j < CWT_OA_SIZE; j++) {
			window[j] = 0;
		}

		// Move data to window buffer and in the data buffer zero out whatever has been moved
		for (int32_t j = 0; j < KERNEL_SIZE; j++) {
			if (step + j < dataLen) {
				window[j * 2] = rral_fixed_to_float(data[step + j]);
				data[step + j] = 0;
			}

			// Add the previous overlap region to the zeroed out segment
			if (step + j + finalOffset >= 0) {
				data[step + j + finalOffset] += overlap[j];
			}

			overlap[j] = 0;
		}

		// Precompute the FFT of the signal segment
		dsps_fft2r_fc32_ansi_(window, CWT_OA_SIZE / 2, fft_table);
		dsps_bit_rev_fc32_ansi(window, CWT_OA_SIZE / 2);

		// Iterate through each scale
		for (uint8_t k = 0; k < CWT_SCALES_SIZE; k++) {

			// Select precalculated kernel from the buffer
			float* kernel = &kernels[CWT_OA_SIZE * k];

			// Perform convolution and return the complex conjugate so the forward DFT can be used as an IDFT 
			for (uint16_t j = 0; (j * 2 + 1) < CWT_OA_SIZE; j++) {
				float rv = window[j * 2] * kernel[j * 2] - window[j * 2 + 1] * kernel[j * 2 + 1];
				float iv = window[j * 2] * kernel[j * 2 + 1] + window[j * 2 + 1] * kernel[j * 2];

				// Complex conjugate of the result
				result[j * 2] = rv;
				result[j * 2 + 1] = -iv;
			}

			// Compute the reverse FFT
			dsps_fft2r_fc32_ansi_(result, CWT_OA_SIZE / 2, fft_table);
			dsps_bit_rev_fc32_ansi(result, CWT_OA_SIZE / 2);

			// Add real results to output
			for (uint16_t j = 0; j < KERNEL_SIZE; j++) {
				if (step + j + finalOffset < dataLen && step + j + finalOffset >= 0) {
					data[step + j + finalOffset] += rral_float_to_fixed(result[j * 2] / (CWT_OA_SIZE / 2));
				}
			}

			// Add real overlap to overlap buffer
			for (uint16_t j = 0; j < KERNEL_SIZE; j++) { // shouldn't this be KERNEL_SIZE-1?
				overlap[j] += rral_float_to_fixed(result[(KERNEL_SIZE + j) * 2] / (CWT_OA_SIZE / 2));
			}
		}
	}

	// Calculate the average value for each point in time
	int32_t d = rral_int_to_fixed(CWT_SCALES_SIZE);
	for (uint16_t i = 0; i < dataLen; i++) {
		if (data[i] != 0) {
			data[i] = rral_div32(data[i], d);
		}
	}

	return RRAL_OK;
}


#endif

