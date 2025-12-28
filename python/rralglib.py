import numpy as np
import math

### The main rralglib equivalent python library ###

### Wrapper function for data validation and simpler algorithm calls
def run_algorithm(data, fs, algorithm="default", args=None):
    
    if len(data) <= 0:
        print("data array cannot be empty")
        return -1, []
    
    if fs <= 0:
        print("sample rate cannot be 0")
        return -1, []
    
    try:
        if algorithm == "default":
            rr, peaks = find_peaks(data=data, fs=fs, args=args)
        elif algorithm == "find_peaks":
            rr, peaks = find_peaks(data=data, fs=fs, args=args)
        elif algorithm == "cwt":
            rr, peaks = cwt_peaks(data=data, fs=fs, args=args)
        elif algorithm == "cwt_oa":
            rr, peaks = cwt_peaks_oa(data=data, fs=fs, args=args)
        elif algorithm == "srmac":
            rr, peaks = srmac(data=data, fs=fs, args=args)
        elif algorithm == "terma":
            rr, peaks = terma(data=data, fs=fs, args=args)
    except Exception as e:
        print("error: "+str(e))
        return -1, []

    return rr, peaks

### Function for finding the approximate respiration window size in samples
def respiration_window(peaks):
    if len(peaks) < 2:
        return 0
    
    first = peaks[0]
    last = peaks[-1]

    ### Copy and reverse input array
    br = [b for b in peaks][::-1]

    ### Calculate the distances between breaths
    for i, b in enumerate(br[:-1]):
        br[i] = br[i] - br[i+1]
    br[-1] = 0
    br = br[::-1]

    ### Calculate the average distance between breaths ### TODO: simply taking an average might not be robust enough; other statistical measures could also be employed
    avg = 0
    for i, b in enumerate(br[1:]):
        avg = (avg * i + b) / (i + 1)

    ### Calculate the respiratory window size
    first = first - avg/2
    last = last + avg/2
    window_size_in_samples = last - first

    return window_size_in_samples

### Function for finding respiratory rate based on breath indices
def find_rr(peaks, fs, actual_window_size):
    peak_count = len(peaks)
    if peak_count < 1:
        return 0
    window_size_in_samples = max(respiration_window(peaks=peaks), actual_window_size*0.75)
    w_time_s = window_size_in_samples/fs
    rr = (60/w_time_s)*peak_count
    return rr

### Function for finding respiratory rate based solely on the average distance between breaths; same concept as above but more efficient
def find_rr_dist(peaks, fs):
    peak_count = len(peaks)
    if peak_count <= 1:
        return 0
    resp_duration = peaks[-1]-peaks[0]
    avg_breath = resp_duration / (peak_count-1)
    rr = (60*fs)/avg_breath
    return rr

### Function for applying a digital biquad iir filter to a signal
def sos_filt(signal, sos):
    
    signal = np.atleast_1d(signal)

    if len(signal) < 3:
        return signal

    ### TODO: initial state needs to be calulated to get rid of the huge transient
    y = np.copy(signal).astype(np.float64)
    for b0, b1, b2, a0, a1, a2 in sos:
        x0, x1, x2 = y[0], y[0], y[0]
        y0, y1, y2 = 0.0, 0.0, 0.0

        for i in range(1, len(y)):
            x0 = y[i]
            y0 = b0*x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2
            y[i] = y0
            x2, x1 = x1, x0
            y2, y1 = y1, y0
    return y

### Function for zero-mean unit variance normalization
def z_norm(signal):
    signal = np.atleast_1d(signal)
    
    if len(signal) == 0:
        return signal
    
    sdev = np.std(signal)
    if sdev != 0:
        return (signal - np.mean(signal)) / sdev
    else:
        return signal - np.mean(signal)
    
### Pearson correlation coefficient
def pearson(x, y, pad=False):
    
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    if len(x) == 0 or len(y) == 0:
        return 0

    ### Pad or truncate arrays if they are not of the same size
    if pad:
        if len(x) != len(y):
            print("x and y arrays are not of the same size. results will be zero-padded to match.")
            max_len = max(len(x), len(y))
            a_ = np.zeros(max_len)
            if len(x) > len(y):
                for i in range(len(y)):
                    a_[i] = y[i]
                y = a_
            else:
                for i in range(len(x)):
                    a_[i] = x[i]
                x = a_
    else:
        if len(x) != len(y):
            print("x and y arrays are not of the same size. results will be truncated to match.")
            min_len = min(len(x), len(y))
            x = np.copy(x[:min_len])
            y = np.copy(y[:min_len])

    ### Calculate the Pearson correlation coefficient
    x_ = np.mean(x)
    y_ = np.mean(y)

    cor = np.sum((x-x_)*(y-y_))
    cor = cor / (math.sqrt(np.sum((x-x_)**2))*math.sqrt(np.sum((y-y_)**2)))

    return cor

### RMSE
def rmse(x, y, pad=True):
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    if len(x) == 0 or len(y) == 0:
        return 0
    
    ### Pad or truncate arrays if they are not of the same size
    if pad:
        if len(x) != len(y):
            print("x and y arrays are not of the same size. results will be zero-padded to match.")
            max_len = max(len(x), len(y))
            a_ = np.zeros(max_len)
            if len(x) > len(y):
                for i in range(len(y)):
                    a_[i] = y[i]
                y = a_
            else:
                for i in range(len(x)):
                    a_[i] = x[i]
                x = a_
    else:
        if len(x) != len(y):
            print("x and y arrays are not of the same size. results will be truncated to match.")
            min_len = min(len(x), len(y))
            x = np.copy(x[:min_len])
            y = np.copy(y[:min_len])

    ### Calculate the RMSE
    rmse = 0
    for i in range(len(x)):
        rmse += abs(x[i]-y[i])**2

    return math.sqrt(rmse / len(x))

### SQI roughly as described in https://ieeexplore.ieee.org/document/6862843
def sqi_full(peaks, raw_signal, max_interval=None):

    peaks = np.atleast_1d(peaks)
    raw_signal = np.atleast_1d(raw_signal)

    if len(peaks) < 2:
        return 0
    
    if len(raw_signal) < 1:
        return 0

    coef = 0
    segments = []

    ### Find the average distance between peaks
    avg_dist = (peaks[-1] - peaks[0]) // (len(peaks))

    ### Check if any distances are larger than max_interval and find min and max distances
    max_dist, min_dist = 0, 0
    if max_interval != None:
        for i in range(len(peaks)-1):
            dist = peaks[i+1] - peaks[i] > max_interval
            if dist > max_interval:
                return 0
            if dist > max_dist:
                max_dist = dist
            elif dist < min_dist:
                min_dist = dist

    ### Save each peak segment
    for peak in peaks:
        if peak - (avg_dist//2) > 0 and peak + (avg_dist//2) < len(raw_signal):
            segments.append(raw_signal[peak-(avg_dist//2):peak+(avg_dist//2)])

    if len(segments) < 1:
        return 0

    segments = np.array(segments)

    ### Calculate the average segment
    avg_segment = np.zeros(len(segments[0]))
    for i in range(len(segments[0])):
        for segment in segments:
            avg_segment[i] += segment[i]
    avg_segment = avg_segment / len(segments)

    ### Calculate the pearson coefficient for each segment
    pearson_coefs = np.zeros(len(segments))
    for i, segment in enumerate(segments):
        pearson_coefs[i] = pearson(segment, avg_segment)

    ### Calculate the average pearson coefficient
    coef = np.mean(pearson_coefs)

    return coef

### SQI naive approach
def sqi_lite(peaks, raw_signal, max_interval=None):

    peaks = np.atleast_1d(peaks)
    raw_signal = np.atleast_1d(raw_signal)

    if len(peaks) < 2:
        return 0
    
    if len(raw_signal) < 1:
        return 0

    coef = 0
    segments = []

    ### Find the average distance between peaks
    avg_dist = (peaks[-1] - peaks[0]) // (len(peaks))

    ### Check if any distances are larger than max_interval and find min and max distances
    max_dist, min_dist = 0, 0
    if max_interval != None:
        for i in range(len(peaks)-1):
            dist = peaks[i+1] - peaks[i] > max_interval
            if dist > max_interval:
                return 0
            if dist > max_dist:
                max_dist = dist
            elif dist < min_dist:
                min_dist = dist

    ### Save each peak segment
    for peak in peaks:
        if peak - (avg_dist//2) > 0 and peak + (avg_dist//2) < len(raw_signal):
            segments.append(raw_signal[peak-(avg_dist//2):peak+(avg_dist//2)])

    if len(segments) < 1:
        return 0

    segments = np.array(segments)

    ### Calculate the pearson coefficient for each neightbouring pair of segments
    pearson_coefs = np.zeros(len(segments)-1)
    for i in range(len(segments)-1):
        pearson_coefs[i] = pearson(segments[i], segments[i+1])

    ### Calculate the average pearson coefficient
    coef = np.mean(pearson_coefs)

    return coef

### Finds all local maxima in data within the range specified by order
def local_maxima_ord(data, order):  

    data = np.atleast_1d(data)

    if len(data) < order*2:
        return []

    maxima = []
    for i in range(order,len(data)-order):
        condition = True
        for j in range(1, order):
            if data[i-j] > data[i]:
                condition = False
                break

            if data[i+j] > data[i]:
                condition = False
                break

        if condition:
            maxima.append(i)

    return maxima

### Finds all local maxima
def local_maxima(data):  

    data = np.atleast_1d(data)

    if len(data) < 3:
        return []
    
    maxima = []
    for i in range(1,len(data)-1):
        condition = True

        if data[i-1] > data[i] or data[i+1] > data[i]:
            condition = False

        if data[i-1] == data[i] and data[i+1] == data[i]:
            condition = False

        if condition:
            maxima.append(i)

    return maxima

### Optimized version of the find_peaks algorithm
def find_peaks(data, fs, window_size=None, prominence=None, heval_ratio=None, width=None, proximity=None, args=None):
    """
    Optimized version of the find_peaks algorithm for peak detection
    """
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []

    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if prominence is None:
            peaks_prom_min = 0.6
        else:
            peaks_prom_min = float(prominence)

        if heval_ratio is None:
            peaks_heval_ratio = 0.8
        else:
            peaks_heval_ratio = float(heval_ratio)

        if width is None:
            peaks_width_min = 0.3 * fs
        else:
            peaks_width_min = float(width) * fs

        if proximity is None:
            peaks_proximity = 1.0
        else:
            peaks_proximity = float(proximity) * fs

    else:
        peaks_prom_min = args[0] if args[0] > -1 else 0.6
        peaks_heval_ratio = args[1] if args[1] > -1 else 0.8
        peaks_width_min = args[2]*fs if args[2] > -1 else 0.3*fs
        peaks_proximity = args[4]*fs if args[4] > -1 else 1.0*fs

    ### Find peaks
    peaks = local_maxima(data)

    ### Return if no peaks found
    if len(peaks) == 0:
        return 0, []
    
    ### Deal with plateaus
    for peak in range(len(peaks)):
        if peak == -1:
            continue
        next = peak + 1
        while next < len(peaks)-1:
            if data[peaks[peak]] == data[peaks[next]]:
                peaks[next] = -1
            else:
                break
            next += 1

    peaks = [peak for peak in peaks if peak != -1]

    ### Remove peaks based on their topographic prominence
    marked = []
    for i, peak in enumerate(peaks):
        neighbour = 0

        ### Find the left prominence
        left_prominence = -1
        for j in range(i-1, -1, -1):
            if j < 0:
                break
            neighbour = peaks[j]
            if data[neighbour] > data[peak]:
                left_prominence = data[peak]-min(data[neighbour:peak])
                break

        if left_prominence == -1:
            left_prominence = data[peak] - min(data[0:peak]) 

        if left_prominence < peaks_prom_min:
            marked.append(i)
            continue
    
        ### Find the right prominence
        right_prominence = -1
        for j in range(i+1, len(peaks)):
            neighbour = peaks[j]
            if data[neighbour] > data[peak]:
                right_prominence = data[peak]-min(data[peak:neighbour])
                break

        if right_prominence == -1:
            right_prominence = data[peak] - min(data[peak:]) 

        if right_prominence < peaks_prom_min:
            marked.append(i)
            continue

        ### Evaluate peak width
        eval_height = data[peak] - peaks_heval_ratio*min(left_prominence, right_prominence)

        width = 0
        for j in range(peak-1, -1, -1):
            if data[j] <= eval_height:
                width += peak - j
                break

        for j in range(peak+1, len(data)):
            if data[j] <= eval_height:
                width += j - peak
                break

        if width < peaks_width_min:
            marked.append(i)
            continue

    ### Remove all marked peaks
    peaks = [peak for i, peak in enumerate(peaks) if i not in marked]
    
    ### Remove peaks that are too close to eachother
    marked = []
    for i in range(len(peaks)):
        if peaks[i] == -1:
            continue
        for j in range(i+1, len(peaks)):
            if peaks[j] == -1:
                continue
            if peaks[j]-peaks[i] >= peaks_proximity:
                break
            else:
                if data[peaks[i]] >= data[peaks[j]]:
                    marked.append(j)
                else:
                    marked.append(i)
                    break
    
    peaks = [peak for i, peak in enumerate(peaks) if i not in marked]

    ### Calculate RR
    # rr = find_rr(peaks, fs, window_size)
    rr = find_rr_dist(peaks, fs)

    return rr, peaks

### Inline SRMAC algorithm ### source: https://arxiv.org/abs/2312.10013
def srmac(data, fs, window_size=None, coef_fast=None, coef_slow=None, coef_cross=None, threshold=None, width=None, margin=None, args=None): 
    """
    SRMAC algorithm for peak detection
    """
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []
    
    if window_size is None:
        window_size = len(data)

    ### Coefficients can be converted from one sampling rate frequency to another using the formula a2 = 1 - (1 - a1)^(f1/f2)
    ### Parameters
    if args is None:
        if coef_fast is None:
            coef_fast = 0.9
        else:
            coef_fast = coef_fast
        
        if coef_slow is None:
            coef_slow = 0.3
        else:
            coef_slow = coef_slow
        
        if coef_cross is None:
            coef_cross = 0.2
        else:
            coef_cross = coef_cross

        if threshold is None:
            th = 0.005
        else:
            th = threshold

        if width is None:
            width = int(0.5 * fs)
        else:
            width = int(width * fs)

        if margin is None:
            margin =  0
        else:
            margin = int(margin * fs)         
    else:
        coef_fast = args[0] if args[0] != -1 else 0.9
        coef_slow = args[1] if args[1] != -1 else 0.3
        coef_cross = args[2] if args[2] != -1 else 0.02
        th = args[3] if args[3] != -1 else 0.005
        width = int(args[4]*fs) if args[4] != -1 else int(0.5*fs)
        margin = int(args[5]*fs) if args[5] != -1 else int(1.0*fs)

    prevdata = data[0]
    prevfast = data[0]
    prevslow = data[0]
    prevcross = 0

    newdata = np.zeros(len(data))

    for i, curr in enumerate(data):
        prevfast = curr * coef_fast + prevfast * (1-coef_fast)
        prevslow = curr * coef_slow + prevslow * (1-coef_slow)
        prevcross = (prevfast-prevslow) * coef_cross + prevcross * (1-coef_cross)
        newdata[i] = prevcross

    peak_count, peaks = zero_crossing(newdata, rawdata=data, width=width, fs=fs, margin=margin, th=th)

    ### Calculate RR
    # rr = find_rr(peaks, fs, window_size)
    rr = find_rr_dist(peaks, fs)

    return rr, peaks

### Corrected version of the TERMA algorithm ### Previous memory and performance optimizations are now quite redundant and need to be rethought
def terma(data, fs, window_size=None, window_event=None, window_cycle=None, b_coef=None, width=None, margin=None, args=None):
    """
    TERMA algorithm for peak detection
    """
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []
    
    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if window_event is None:
            w1 = int(1 * fs)
        else:
            w1 = int(window_event * fs)

        if window_cycle is None:
            w2 = int(3 * fs)
        else:
            w2 = int(window_cycle * fs)

        if b_coef is None:
            b = 0.5
        else:
            b = b_coef

        if width is None:
            width = int(0.7 * fs)
        else:
            width = int(width * fs)

        if margin is None:
            margin =  0
        else:
            margin = int(margin * fs) 
    else:
        w1 = int(args[0]*fs) if args[0] != -1 else fs
        w2 = int(args[1]*fs) if args[1] != -1 else 3*fs
        b = args[2] if args[2] != -1 else 0.5
        width = int(args[3]*fs) if args[3] != -1 else int(0.7*fs)
        margin = int(args[4]*fs) if args[4] != -1 else int(1.0*fs)

    ### Clone input buffer so the input data remains unchanged
    data_ = np.array(data)
    data_event = np.array(data)
    data_cycle = np.array(data)

    circBuf_event = np.zeros(w1)
    circBuf_cycle = np.zeros(w2)
    wptr_ev = 0
    wptr_cy = 0

    ev_sum = 0
    ev_oldest = 0
    ev_newest = 0
    
    cy_sum = 0
    cy_oldest = 0
    cy_newest = 0

    # Check if circular buffers are larger than input data
    if len(circBuf_cycle) > len(data_) or len(circBuf_event) > len(data):
        return 0, []

    # Starting state of event circular buffer
    for i in range(len(circBuf_event)):
        if i < int(w1/2):
            circBuf_event[i] = 0
        else:
            circBuf_event[i] = data_[i-int(w1/2)]
            ev_sum += circBuf_event[i]
        
    # Starting state of cycle circular buffer
    for i in range(len(circBuf_cycle)):
        if i < int(w2/2):
            circBuf_cycle[i] = 0
        else:
            circBuf_cycle[i] = data_[i-int(w2/2)]
            cy_sum += circBuf_cycle[i]
        
    z = np.mean(data_)

    # Iterate over all elements in data
    for i in range(0, len(data_)):
        event = ev_sum/w1
        cycle = cy_sum/w2

        # Calculate new data value
        data_event[i] = event
        data_cycle[i] = cycle + b * z
        data_[i] = data_event[i] - data_cycle[i]

        # Update event buffer data
        if i+int(w1/2) >= len(data_):
            ev_oldest = circBuf_event[wptr_ev]
            circBuf_event[wptr_ev] = 0
            ev_newest = 0
        else:    
            ev_oldest = circBuf_event[wptr_ev]
            circBuf_event[wptr_ev] = data_[i+int(w1/2)]
            ev_newest = circBuf_event[wptr_ev]

        # Increment event write pointer
        wptr_ev = (wptr_ev + 1) % w1

        # Update cycle buffer data
        if i+int(w2/2) >= len(data_):
            cy_oldest = circBuf_cycle[wptr_cy]
            circBuf_cycle[wptr_cy] = 0
            cy_newest = 0
        else:    
            cy_oldest = circBuf_cycle[wptr_cy]
            circBuf_cycle[wptr_cy] = data_[i+int(w2/2)]
            cy_newest = circBuf_cycle[wptr_cy]

        # Increment cycle write pointer
        wptr_cy = (wptr_cy + 1) % w2

        # Recalculate sums
        ev_sum = ev_sum - ev_oldest + ev_newest
        cy_sum = cy_sum - cy_oldest + cy_newest

    # fig, ax = plt.subplots(1,1)
    # # ax.plot(data_event)
    # # ax.plot(data_cycle)
    # ax.plot(data_event-data_cycle, c='r')
    # plt.show()

    peak_count, peaks = zero_crossing(data=data_, rawdata=data, width=w1, fs=fs, margin=margin, th=0.0)

    ### Calculate RR
    # rr = find_rr(peaks, fs, window_size)
    rr = find_rr_dist(peaks, fs)

    return rr, peaks

### Count-orig peak detection method described in https://link.springer.com/article/10.1007/s10439-007-9428-1 ### https://github.com/peterhcharlton/RRest/blob/master/RRest_v3.0/Algorithms/estimate_rr/CtO.m
def count_orig(data, fs, window_size=None, max_troughs=None, percentile=None, th_coef=None, args=None):

    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []
    
    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if max_troughs is None:
            max_troughs = 1
        else:
            max_troughs = int(max_troughs)

        if percentile is None:
            percentile = 75
        else:
            percentile = int(percentile)

        if th_coef is None:
            th_coef = 0.2
        else:
            th_coef = float(th_coef)
    else:
        max_troughs = int(args[0]) if args[0] != -1 else 1
        percentile = int(args[1]) if args[1] != -1 else 75
        th_coef = float(args[2]) if args[2] != -1 else 0.2

    ### Find local maxima
    peaks = local_maxima(data)
    troughs = local_maxima(-data)

    if len(peaks) == 0 or len(troughs) == 0:
        return 0, []

    ### Define threshold
    q3 = np.percentile(data[peaks], percentile)
    th = th_coef * q3
    
    ### Find relevant peaks and troughs
    peaks = [peak for peak in peaks if data[peak] > th]
    troughs = [trough for trough in troughs if data[trough] < 0]

    ### Find valid respiratory cycles
    cycles = []
    cycle_durations = []
    for i in range(len(peaks)-1):
        cycle_troughs = 0
        ### A cycle is only considered valid if there is one trough between two cycles - this can be changed with the max_troughs parameter
        for trough in troughs:
            if trough > peaks[i] and trough < peaks[i+1]:
                cycle_troughs += 1
        if cycle_troughs <= max_troughs and cycle_troughs > 0:
            cycles.append(peaks[i])
            cycle_durations.append(peaks[i+1]-peaks[i])

    ### Find RR
    mean_duration = np.mean(cycle_durations)
    if mean_duration == 0:
        rr = 0
    else:
        rr = (60*fs)/mean_duration

    return rr, cycles

### Count-adv peak detection method described in https://link.springer.com/article/10.1007/s10439-007-9428-1 ### https://github.com/peterhcharlton/RRest/blob/master/RRest_v3.0/Algorithms/estimate_rr/CtA.m
def count_adv(data, fs, window_size=None, percentile=None, th_coef=None, args=None): 
    
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []

    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if percentile is None:
            percentile = 75
        else:
            percentile = int(percentile)

        if th_coef is None:
            th_coef = 0.8
        else:
            th_coef = float(th_coef)
    else:
        percentile = int(args[1]) if args[1] != -1 else 75
        th_coef = float(args[2]) if args[2] != -1 else 0.8

    ### Find local maxima
    peaks = local_maxima(data)
    troughs = local_maxima(-data)

    ### Define threshold
    extrema = peaks[:]
    for trough in troughs:
        extrema.append(trough)
    if len(extrema) == 0:
        return 0, []
    extrema = np.sort(extrema)
    amps = data[extrema]
    amp_diffs = abs(amps[1:-1]-amps[2:])
    q3 = np.percentile(amp_diffs, percentile)
    th = th_coef * q3

    ### Eliminate pairs of extrema if the difference in amplitude between them is smaller than the threshold
    eliminating = 1
    while eliminating:
        if len(extrema) < 3:
            eliminating = 0
            continue
        amps = data[extrema]
        amp_diffs = amp_diffs = abs(amps[:-1]-amps[1:])

        min_amp_diff_idx = np.argmin(amp_diffs)
        if amp_diffs[min_amp_diff_idx] > th:
            eliminating = 0
        else:
            min_amp_pair_idx = min_amp_diff_idx + 1
            extrema = [e for i, e in enumerate(extrema) if i != min_amp_diff_idx and i != min_amp_pair_idx]

    if len(extrema) < 3:
        return 0, []

    ### Truncate to start and end at a maximum
    if data[extrema[0]] < data[extrema[1]]:
        extrema = extrema[1:]
    if data[extrema[-1]] < data[extrema[-2]]:
        extrema = extrema[:-1]

    ### Find RR
    no_breaths = (len(extrema)-1)/2
    breathing_duration = extrema[-1]-extrema[0]
    avg_breath_duration = breathing_duration / no_breaths
    rr = (60*fs)/avg_breath_duration

    return rr, extrema[::2]

### Method for finding peaks in zero-crossing data
def zero_crossing(data, width, fs=None, th=0, rawdata=None, margin=None):
    
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []

    if len(data) < margin*2:
        return 0, []

    if rawdata is None:
        rawdata = data

    if margin is None:
        margin = 0

    if width < 1:
        width = 1

    positive, peak_count, delta = 0, 0, 0
    maximum = th
    peaks = []

    for i in range(margin, len(data)-margin):
        if data[i] > th:
            positive += 1
            delta += 1
            if rawdata[i] > maximum:
                maximum = rawdata[i]
                delta = 0
        else:
            delta += 1
            if positive >= width:
                peaks.append(i-delta)
                peak_count += 1
            positive, delta = 0, 0
            maximum = th

    return peak_count, peaks

### Mexh wavelet 
def mexh(t): ### https://pywavelets.readthedocs.io/en/latest/ref/cwt.html#mexican-hat-wavelet
    return ((2/3**(1/2))*((1/math.pi)**(1/4)))*math.exp(-t**2/2)*(1-t**2)

### Second order Gaussian wavelet
def gaus2(t): ### https://ieeexplore.ieee.org/document/9941083
    return ((2/3**(1/2))*((2/math.pi)**(1/4)))*math.exp(-t**2)*(1-2*t**2)

### Wrapper function for scaled wavelets ### https://ccrma.stanford.edu/~jos/sasp/Continuous_Wavelet_Transform.html
def wavelet(t, tau, scale, wavelet="gaus2"): 
    if scale == 0:
        print("scale cannot be 0")
        return 0
    if wavelet == "gaus2":
        return (1/scale**(1/2))*gaus2((t-tau)/scale)
    elif wavelet == "mexh":
        return (1/scale**(1/2))*mexh((t-tau)/scale)
    else:
        print(wavelet + " wavelet is not a valid option.")
        return 0

### Index bit-reverse sort algorithm ### Required as part of the FFT procedure in order for results to be returned in natural order
def bit_reverse_sort(x):  ### Pseudocode source: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Data_reordering,_bit_reversal,_and_in-place_algorithms

    x = np.atleast_1d(x)

    N = len(x)

    if math.log2(N) % 1 != 0.0:
        print("error: cannot perform bit-reverse-sort on arrays not of size log2")
        return x

    bits = int(math.log2(N))

    y = np.zeros(N, dtype=complex)

    for k in range(N):
        k_b = bin(k)[2:].zfill(bits)
        k_b = int(k_b[::-1], 2)
        y[k_b] = x[k]

    return y

### Iterative radix-2 DIT FFT algorithm
def fft(x):
    N = len(x)

    if N == 1:
        print("Input data is of len=1")
        return x

    if math.log2(N) % 1 != 0.0:
        print("Input data not of len=2**p")
        return x

    ### Cast array to complex numbers
    x = np.array(x, dtype=complex)

    ### Bit-reverse-copy procedure so results are in natural order
    x = bit_reverse_sort(x=x)

    for s in range(1,int(math.log2(N))+1): ### Pseudocode source: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Data_reordering,_bit_reversal,_and_in-place_algorithms
        m = 2**s
        wm = math.e**((-2j*math.pi)/m)
        for k in range(0, N, m):
            w = 1
            for j in range(m//2):
                t = w*x[k + j + m//2]
                u = x[k + j]
                x[k + j] = (u + t)
                x[k + j + m//2] = (u - t)
                w = w*wm
    return x

### Inverse FFT
def ifft(x): 
    N = len(x)

    if N == 1:
        print("Input data is of len=1")
        return x

    if math.log2(N) % 1 != 0.0:
        print("Input data not of len=2**p")
        return x
    
    ### Cast array to complex numbers
    x = np.array(x, dtype=complex)

    ### Bit-reverse-copy procedure so results are in natural order
    x = bit_reverse_sort(x=x)

    for s in range(1,int(math.log2(N))+1):
        m = 2**s
        wm = math.e**(2j*math.pi/m)
        for k in range(0, N, m):
            w = 1
            for j in range(m//2):
                t = w*x[k + j + m//2]
                u = x[k + j]
                x[k + j] = (u + t)
                x[k + j + m//2] = (u - t)
                w = w*wm
    return x/N ### Output is scaled by N

### Pad array with zeroes up to len(x) == n
def padn(x, n):
    if n <= len(x):
        return x
    xp = np.zeros(int(n))
    for i in range(len(x)):
        xp[i] = x[i]
    return xp

### Custom CWT peaks implementation ### Reference: https://pywavelets.readthedocs.io/en/latest/ref/cwt.html#
def cwt_peaks(data, fs, window_size=None, resolution=None, threshold=None, width=None, margin=None, min_freq=None, max_freq=None, kernel_size=None, args=None):
    """
    CWT algorithm for peak detection
    """
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []

    if fs <= 0:
        return 0, []

    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if resolution is None:
            res = 15
        else:
            res = int(resolution)
        
        if threshold is None:
            th = 0.1
        else:
            th = float(threshold)

        if width is None:
            width = int(0.5 * fs)
        else:
            width = int(width * fs)

        if margin is None:
            margin =  0
        else:
            margin = int(margin * fs)

        if min_freq is None:
            f_min = 0.02
        else:
            f_min = min_freq

        if max_freq is None:
            f_max = 0.73
        else:
            f_max = max_freq

        if kernel_size is None:
            wavelet_length = 2
        else:
            wavelet_length = kernel_size
    else:
        res = int(args[0]) if args[0] != -1 else 15
        th = args[1] if args[1] != -1 else 0.1
        width = int(args[2]*fs) if args[2] != -1 else int(0.5*fs)
        margin = int(args[3]*fs) if args[3] != -1 else int(1.0*fs)
        f_min = args[4] if args[4] != -1 else 0.02
        f_max = args[5] if args[5] != -1 else 0.73
        wavelet_length = kernel_size

    res = 5
    centfreq = 0.25 ### https://scispace.com/pdf/qrs-complex-detection-using-combination-of-mexican-hat-30mk07j6d6.pdf

    ### Calculate appropriate CWT scales
    scales = centfreq/(np.linspace(start=f_min, stop=f_max, num=res)/fs)
    
    ### Base array for the wavelet kernels
    w_t = np.arange(start=-(wavelet_length/2)*fs, stop=(wavelet_length/2)*fs)

    ### Pad the input data with zeroes up to len(data) == nearest power of 2
    ### TODO: formalize this process
    data_len = len(data)
    kern_len = len(w_t)
    size = data_len+kern_len-1

    p = math.log2(size)
    if p % 1 != 0.0:
        p = int(p)+1

    size = 2**p

    data = padn(data, size)
    
    ### For each scale create a discretized and scaled version of the wavelet and perform convolution on the signal
    ### TODO: as long as the algorithm parameters are unchanged, wavelets(and even their FFTs) do not need be recalculated and could be stored and reused for a slight performance gain
    out = []
    for i, scale in enumerate(scales):
        kernel = []
        for i in range(len(w_t)):
            kernel.append(wavelet(w_t[i], 0, scale, "gaus2"))

        ### Pad the kernel with zeroes to match len(data)
        kernel = padn(kernel, size)

        ### TODO: Overlap-add convolution (is not as good for small fft sizes)
        out.append(ifft(fft(data)*fft(kernel))[int(kern_len/2):data_len+int(kern_len/2)+1]) #

    ### Average across scales
    data_avg = []
    for i in range(len(out[0])):
        val = 0
        for scale in range(res):
            # val += out[scale][i]/scales[scale] ### normalize power (source?)
            val += out[scale][i]

        val = val/res
        data_avg.append(val)

    # return data_avg

    ### Find breaths
    peak_count, peaks = zero_crossing(data_avg, width=width, fs=fs, margin=margin, th=th)

    ### Calculate RR
    # rr = find_rr(peaks, fs, window_size)
    rr = find_rr_dist(peaks, fs)

    return rr, peaks

### Overlap-add version of the CWT peaks algorithm
def cwt_peaks_oa(data, fs, window_size=None, resolution=None, threshold=None, width=None, margin=None, min_freq=None, max_freq=None, kernel_size=None, args=None):
    """
    Overlap-add CWT algorithm for peak detection
    """
    data = np.atleast_1d(data)

    if len(data) < 1:
        return 0, []
    
    if fs <= 0:
        return 0, []

    if window_size is None:
        window_size = len(data)

    ### Parameters
    if args is None:
        if resolution is None:
            res = 15
        else:
            res = int(resolution)
        
        if threshold is None:
            th = 0.0
        else:
            th = float(threshold)

        if width is None:
            width = int(0.5 * fs)
        else:
            width = int(width * fs)

        if margin is None:
            margin =  0
        else:
            margin = int(margin * fs)

        if min_freq is None:
            f_min = 0.02
        else:
            f_min = min_freq

        if max_freq is None:
            f_max = 0.73
        else:
            f_max = max_freq

        if kernel_size is None:
            wavelet_length = 2
        else:
            wavelet_length = kernel_size
    else:
        res = int(args[0]) if args[0] != -1 else 15
        th = args[1] if args[1] != -1 else 0.0
        width = int(args[2]*fs) if args[2] != -1 else int(0.5*fs)
        margin = int(args[3]*fs) if args[3] != -1 else int(1.0*fs)
        f_min = args[4] if args[4] != -1 else 0.02
        f_max = args[5] if args[5] != -1 else 0.73
        if kernel_size is None:
            wavelet_length = 2
        else:
            wavelet_length = kernel_size

    centfreq = 0.25 ### https://scispace.com/pdf/qrs-complex-detection-using-combination-of-mexican-hat-30mk07j6d6.pdf

    ### Calculate appropriate CWT scales
    scales = centfreq/(np.linspace(start=f_min, stop=f_max, num=res)/fs)
    
    ### Base array for the wavelet kernels
    wavelet_t = np.arange(start=-(wavelet_length/2)*fs, stop=(wavelet_length/2)*fs)

    ### Calculate the size of overlap-add FFT windows
    signal_len = len(data)
    kernel_len = wavelet_length * fs
    conv_size = (kernel_len * 2) - 1

    p = math.log2(conv_size)
    if p % 1 != 0.0:
        p = int(p)+1

    conv_size = 2**p

    ### Divide signal into equal parts
    window_start = 0
    signal_windows = []
    while window_start < signal_len:
        temp = np.zeros(kernel_len)
        for i in range(kernel_len):
            if window_start + i < signal_len:
                temp[i] = data[window_start + i]
            else:
                temp[i] = 0

        ### Pad the window with zeroes for FFT convolution
        signal_windows.append(padn(temp, conv_size))
            
        ### Move window forward
        window_start += kernel_len

    signal_windows = np.array(signal_windows)

    ### For each scale create a scaled version of the wavelet
    kernels = []
    for i, scale in enumerate(scales):
        kernel = np.zeros(kernel_len)
        for i in range(kernel_len):
            kernel[i] = wavelet(wavelet_t[i], 0, scale, "gaus2") # * (1/math.sqrt(scale))

        ### Pad the kernel with zeroes for FFT convolution
        kernels.append(padn(kernel, conv_size))

    kernels = np.array(kernels)

    ### For each kernel perform overlap-add FFT convolution over all signal windows
    scalogram = []
    for kernel in kernels:
        convolved_signal = np.zeros(kernel_len*len(signal_windows) + int(1.5 * kernel_len) + 1)
        for i, window in enumerate(signal_windows):
            convolved_window = ifft(fft(window)*fft(kernel))[:kernel_len*2]
            for j in range(len(convolved_window)):
                convolved_signal[i * kernel_len + j] += convolved_window[j].real

        ### Slicing performed to remove phase shift
        scalogram.append(convolved_signal[int(kernel_len/2):int(-kernel_len/2)-fs-1])

    ### Average across scales
    data_avg = np.zeros(signal_len)
    for i in range(signal_len):
        val = 0
        for scale in range(res):
            val += scalogram[scale][i]

        val = val/res
        data_avg[i] = val

    ### Find breaths
    peak_count, peaks = zero_crossing(data_avg, width=width, fs=fs, margin=margin, th=th)

    ### Calculate RR
    # rr = find_rr(peaks, fs, window_size)
    rr = find_rr_dist(peaks, fs)

    return rr, peaks
