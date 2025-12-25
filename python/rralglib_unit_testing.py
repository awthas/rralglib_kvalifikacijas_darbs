import numpy as np
import unittest
import rr_algorithms

### Prototype unit tests for the Python rralglib module

class TestSRMAC(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rr_algorithms.srmac_inline(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rr_algorithms.srmac_inline(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rr_algorithms.srmac_inline(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rr_algorithms.srmac_inline(data=signal, fs=fs, args=[-1,-1,-1,0,0.1,0])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestTERMA(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rr_algorithms.terma_corrected(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rr_algorithms.terma_corrected(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rr_algorithms.terma_corrected(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rr_algorithms.terma_corrected(data=signal, fs=fs, args=[0.2,0.6,-1,0.1,0])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestFindPeaks(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rr_algorithms.find_peaks_fast(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rr_algorithms.find_peaks_fast(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rr_algorithms.find_peaks_fast(data=signal, fs=1)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rr_algorithms.find_peaks_fast(data=signal, fs=fs, args=[-1,-1,-1,-1,-1,-1])

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestCWTPeaks(unittest.TestCase):

    def test_empty(self):
        signal = []

        rr, peaks = rr_algorithms.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = 1

        rr, peaks = rr_algorithms.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        rr, peaks = rr_algorithms.cwt_peaks_oa(data=signal, fs=1, threshold=0.0)

        self.assertEqual(rr, 0)
        self.assertEqual(peaks, [])

    def test_sine(self):
        fs = 64
        signal = np.sin(np.linspace(0,4*np.pi,fs*4))

        rr, peaks = rr_algorithms.cwt_peaks_oa(data=signal, fs=fs, resolution=5, threshold=0.0, width=0.1, margin=0, min_freq=0.5, max_freq=4)

        self.assertGreater(rr, 0)
        self.assertGreater(len(peaks), 0)

class TestZeroCrossing(unittest.TestCase):

    def test_empty(self):
        signal = []

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_one(self):
        signal = [1]

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_zeroes(self):
        signal = np.zeros(100)

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 0)
        self.assertEqual(peaks, [])

    def test_alternating(self):
        signal = np.array([0,1,0,1,0,1,0,1,0])

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertEqual(count, 4)
        self.assertEqual(peaks, [1,3,5,7])

    def test_alternating_threshold(self):
        signal = np.array([0,0.5,0,0.7,0,0.9,0,1.1,0])

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.8)

        self.assertEqual(count, 2)
        self.assertEqual(peaks, [5,7])

    def test_sine(self):
        signal = np.sin(np.linspace(0,4*np.pi,256))

        count, peaks = rr_algorithms.zero_crossing(signal, width=1, fs=0, rawdata=None, margin=0, th=0.0)

        self.assertGreater(count, 0)
        self.assertGreater(len(peaks), 0)

class TestFFT(unittest.TestCase):

    def test_zeroes(self):
        signal = np.zeros(256)
        signal = rr_algorithms.fft(signal)

        self.assertListEqual(signal.tolist(), np.zeros(256).tolist())

    def test_sine_single(self):
        signal = np.sin(np.linspace(0,16*np.pi,256)) ### 1Hz sine wave at a sampling rate of 32Hz
        signal = rr_algorithms.fft(signal) ### n * 32 / 256 = 1Hz, therefore peak is expected at index 8

        self.assertGreater(signal[8], signal[7])
        self.assertGreater(signal[8], signal[9])

    def test_sine_composite(self):
        signal = np.sin(np.linspace(0,16*np.pi,256))+np.sin(np.linspace(0,64*np.pi,256)) ### 1Hz sine wave and 4Hz sine wave at a sampling rate of 32Hz
        signal = rr_algorithms.fft(signal) ### n * 32 / 256 = 1Hz/4Hz, therefore peaks expected at index 8 and 32

        self.assertGreater(signal[8], signal[7])
        self.assertGreater(signal[8], signal[9])
        self.assertGreater(signal[32], signal[31])
        self.assertGreater(signal[32], signal[33])

class TestIFFT(unittest.TestCase):

    def test_zeroes(self):
        signal = np.zeros(256)
        fftsignal = rr_algorithms.fft(signal)
        fftsignal = rr_algorithms.ifft(fftsignal)

        self.assertListEqual(signal.tolist(), fftsignal.real.tolist())

    def test_sine(self):
        signal = np.sin(np.linspace(0,16*np.pi,256)) ### 1Hz sine wave at a sampling rate of 32Hz
        fftsignal = rr_algorithms.fft(signal)
        fftsignal = rr_algorithms.ifft(fftsignal).real

        signal = signal.tolist()
        for i in range(len(signal)):
            signal[i] = round(signal[i], 8)

        fftsignal = fftsignal.tolist()
        for i in range(len(fftsignal)):
            fftsignal[i] = round(fftsignal[i], 8)
        
        self.assertListEqual(signal, fftsignal)

# sosfilt

# pearson

# rmse

# sqi

# local maxima

# count_orig

# count_adv

# bit rev

if __name__ == '__main__':
    unittest.main()


