import numpy as np
import unittest
import rralglib

### Unit tests for the Python rralglib module

class TestSRMAC(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rralglib.srmac(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.srmac(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.srmac(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.srmac(data=signal, fs=fs, args=[-1,-1,-1,0,0.1,0])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestTERMA(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rralglib.terma(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.terma(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.terma(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.terma(data=signal, fs=fs, args=[0.2,0.6,-1,0.1,0])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestFindPeaks(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rralglib.find_peaks(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.find_peaks(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.find_peaks(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.find_peaks(data=signal, fs=fs, args=[-1,-1,-1,-1,-1,-1])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestCWT(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rralglib.cwt_peaks(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.cwt_peaks(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.cwt_peaks(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.cwt_peaks(data=signal, fs=fs, resolution=5, threshold=0.0, width=0.1, margin=0, min_freq=0.5, max_freq=4)

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestCWTOA(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rralglib.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.cwt_peaks_oa(data=signal, fs=fs, resolution=5, threshold=0.0, width=0.1, margin=0, min_freq=0.5, max_freq=4)

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestZeroCrossing(unittest.TestCase):

    def test_empty(self):
        signal = []

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = [1]

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_alternating(self):
        signal = np.array([0,1,0,1,0,1,0,1,0])

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 4)
        self.assertEqual(peaks, [1,3,5,7])

    def test_alternating_threshold(self):
        signal = np.array([0,0.5,0,0.7,0,0.9,0,1.1,0])

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.8)

        self.assertEqual(count, 2)
        self.assertEqual(peaks, [5,7])

    def test_sine(self):
        signal = np.sin(np.linspace(0,4*np.pi,256))

        count, peaks = rralglib.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertGreater(count, 0)
        self.assertGreater(len(peaks), 0)

class TestFFT(unittest.TestCase):

    def test_zeroes(self):
        signal = np.zeros(256)
        signal = rralglib.fft(signal)

        self.assertListEqual(signal.tolist(), np.zeros(256).tolist())

    def test_sine_single(self):
        signal = np.sin(np.linspace(0,16*np.pi,256)) ### 1Hz sine wave at a sampling rate of 32Hz
        signal = rralglib.fft(signal) ### n * 32 / 256 = 1Hz, therefore peak is expected at index 8

        self.assertGreater(signal[8], signal[7])
        self.assertGreater(signal[8], signal[9])

    def test_sine_composite(self):
        signal = np.sin(np.linspace(0,16*np.pi,256))+np.sin(np.linspace(0,64*np.pi,256)) ### 1Hz sine wave and 4Hz sine wave at a sampling rate of 32Hz
        signal = rralglib.fft(signal) ### n * 32 / 256 = 1Hz/4Hz, therefore peaks expected at index 8 and 32

        self.assertGreater(signal[8], signal[7])
        self.assertGreater(signal[8], signal[9])
        self.assertGreater(signal[32], signal[31])
        self.assertGreater(signal[32], signal[33])

class TestIFFT(unittest.TestCase):

    def test_zeroes(self):
        signal = np.zeros(256)
        fftsignal = rralglib.fft(signal)
        fftsignal = rralglib.ifft(fftsignal)

        self.assertListEqual(signal.tolist(), fftsignal.real.tolist())

    def test_sine(self):
        signal = np.sin(np.linspace(0,16*np.pi,256)) ### 1Hz sine wave at a sampling rate of 32Hz
        fftsignal = rralglib.fft(signal)
        fftsignal = rralglib.ifft(fftsignal).real

        signal = signal.tolist()
        for i in range(len(signal)):
            signal[i] = round(signal[i], 8)

        fftsignal = fftsignal.tolist()
        for i in range(len(fftsignal)):
            fftsignal[i] = round(fftsignal[i], 8)
        
        self.assertListEqual(signal, fftsignal)

class TestWavelet(unittest.TestCase): ### All assert values calculated using a calculator and the appropriate wavelet definitions
    def test_gaus2_zero(self):
        t = 0

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="gaus2")

        self.assertEqual(1.0314, round(wavelet, 4))

    def test_gaus2_one(self):
        t = 1

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="gaus2")

        self.assertEqual(-0.3794, round(wavelet, 4))

    def test_gaus2_pi(self):
        t = np.pi

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="gaus2")

        self.assertEqual(-0.0010, round(wavelet, 4))

    def test_mexh_zero(self):
        t = 0

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="mexh")

        self.assertEqual(0.8673, round(wavelet, 4))

    def test_mexh_one(self):
        t = 1

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="mexh")

        self.assertEqual(0.0, round(wavelet, 4))

    def test_mexh_pi(self):
        t = np.pi

        wavelet = rralglib.wavelet(t, tau=0, scale=1, wavelet="mexh")

        self.assertEqual(-0.0553, round(wavelet, 4))

class TestSosfilt(unittest.TestCase): 
    def test_lowpass_sine(self): ### Sine wave at a frequency of 1Hz, lowpass filtered with a 3rd order Butterworth filter at 0.5Hz; the signal after the transient should be close to zero
        fs = 4
        signal = np.sin(np.linspace(0,256*np.pi,fs*64))
        sos = [[0.03168934,0.06337869,0.03168934,1.,-0.41421356,0.],[1.,1.,0.,1.,-1.0448155,0.47759225]]

        signal = rralglib.sos_filt(signal, sos)

        self.assertListEqual(np.zeros(64).tolist(), np.round(signal[-65:-1], 5).tolist())

class TestPearson(unittest.TestCase):

    def test_twonumbers_equal(self):
        x = 1
        y = 1

        r = rralglib.pearson(x, y)

        self.assertAlmostEqual(1.0, r)

    def test_twoarrays_equal(self):
        x = np.arange(128)
        y = np.arange(128)

        r = rralglib.pearson(x, y)

        self.assertAlmostEqual(1.0, r)

    def test_twoarrays_opposite(self):
        x = np.arange(128) - 64
        y = np.copy(x) * -1

        r = rralglib.pearson(x, y)

        self.assertAlmostEqual(-1.0, r)

class TestSQIFull(unittest.TestCase):
    def test_onepeak(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([20])

        sqi = rralglib.sqi_full(peaks, signal)

        self.assertAlmostEqual(0.0, sqi)

    def test_fourpeaks(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([20, 28, 36, 52])

        sqi = rralglib.sqi_full(peaks, signal)

        self.assertAlmostEqual(1.0, sqi)

    def test_twopeaks_unequal(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([24, 48])

        signal[25] = 0
        signal[23] = 0

        sqi = rralglib.sqi_full(peaks, signal)

        self.assertGreater(1.0, sqi)

class TestSQILite(unittest.TestCase):
    def test_onepeak(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([20])

        sqi = rralglib.sqi_lite(peaks, signal)

        self.assertAlmostEqual(0.0, sqi)

    def test_fourpeaks(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([20, 28, 36, 52])

        sqi = rralglib.sqi_lite(peaks, signal)

        self.assertAlmostEqual(1.0, sqi)

    def test_twopeaks_unequal(self):
        pattern = np.array([1,2,3,4,5,4,3,2])
        signal = np.zeros(128)

        for i in range(0,len(signal),len(pattern)):
            signal[i:i+len(pattern)] = pattern

        peaks = np.array([24, 48])

        signal[25] = 0
        signal[23] = 0

        sqi = rralglib.sqi_lite(peaks, signal)

        self.assertGreater(1.0, sqi)

class TestCountOrig(unittest.TestCase):
    def test_empty(self):
        signal = []

        rr, peaks = rralglib.count_orig(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.count_orig(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.count_orig(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.count_orig(data=signal, fs=fs)

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestCountAdv(unittest.TestCase):
    def test_empty(self):
        signal = []

        rr, peaks = rralglib.count_adv(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rralglib.count_adv(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rralglib.count_adv(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rralglib.count_adv(data=signal, fs=fs)

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestBitrev(unittest.TestCase): ### Examples given in https://en.wikipedia.org/wiki/Bit-reversal_permutation used as test cases
    def test_zero(self):
        x = np.array([0])

        x = rralglib.bit_reverse_sort(x)

        self.assertListEqual([0],x.tolist())

    def test_range_two(self):
        x = np.arange(2)

        x = rralglib.bit_reverse_sort(x)

        self.assertListEqual([0,1],x.tolist())

    def test_range_four(self):
        x = np.arange(4)

        x = rralglib.bit_reverse_sort(x)

        self.assertListEqual([0,2,1,3],x.tolist())

    def test_range_sixteen(self):
        x = np.arange(16)

        x = rralglib.bit_reverse_sort(x)

        self.assertListEqual([0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15],x.tolist())

if __name__ == '__main__':
    unittest.main()


