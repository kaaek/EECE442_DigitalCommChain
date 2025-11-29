# Digital Communication Chain

**Authors:** Khalil El Kaaki, Mouhammad Kandakji  
**Course:** EECE 442: Communication Systems  
**Institution:** American University of Beirut

---

## Overview

This project implements a comprehensive digital communication system spanning:
- **Phase 1:** Signal sampling, Fourier series analysis, quantization, and lossless source coding
- **Phase 2:** BPSK/QPSK modulation/demodulation, channel simulation (AWGN, ISI), and MLSE equalization with BER analysis

---

## Repository Structure

### Main Entry Point
- **`main.m`** - Master script controlling all experiment modules (Phase 1 & Phase 2)
- **`constants.m`** - Global configuration: signal duration, frequencies, MSE targets, quantizer levels

### Phase 1: Analog Signal Processing

#### Basic Signal Generation
- **`Analog/AMWave.m`** - Generates AM-modulated signal
- **`Analog/sample.m`** - Samples continuous signal at specified frequency Fs
- **`Analog/reconstruct.m`** - Reconstructs signal from samples using sinc interpolation
- **`Analog/sinc.m`** - Sinc function: sin(π*x)/(π*x), with limit as x→0 = 1

#### Fourier Series Analysis
- **`Analog/Fourier/fcoef.m`** - Computes j-th Fourier series coefficient via numerical integration
- **`Analog/Fourier/fs.m`** - Computes Fourier series approximation with N coefficients
- **`Analog/Fourier/fsConstN.m`** - Analyzes approximation with fixed N, varying period T
- **`Analog/Fourier/fsConstT.m`** - Analyzes convergence with fixed period T, varying N

### Phase 2: Modulation/Demodulation

#### BPSK (Binary Phase Shift Keying)
- **`Analog/PSK/bpsk_mod.m`** - Modulates bit (0,1) to constellation: 0→-A, 1→+A
- **`Analog/PSK/bpsk_demod.m`** - Demodulates symbol via threshold detection at 0

#### QPSK (Quadrature Phase Shift Keying)
- **`Analog/PSK/qpsk_mod.m`** - Modulates 2-bit string to QPSK constellation (counter-clockwise)
- **`Analog/PSK/qpsk_demod.m`** - Demodulates symbols via quadrant-based detection

#### Utility
- **`Analog/PSK/random_bitstream.m`** - Generates random bit sequence (default 10e6 bits)

### Phase 1: Digital Signal Processing

#### Quantization
- **`Digital/Quantizer/quantize.m`** - Maps samples to quantization levels via thresholds
- **`Digital/Quantizer/partition.m`** - Computes optimal region centroids

**Uniform Quantizer:**
- **`Digital/Quantizer/Uniform Quantizer/uniformQuan.m`** - M equally-spaced quantization levels

**Two-Level Quantizer:**
- **`Digital/Quantizer/Uniform Quantizer/twoLvlQuan.m`** - Binary quantization at mean threshold

**Lloyd-Max Quantizer (Optimal):**
- **`Digital/Quantizer/Lloyd-Max Quantizer/lloydMax.m`** - Iterative optimization algorithm with target MSE convergence
- **`Digital/Quantizer/Lloyd-Max Quantizer/lloydMaxInit.m`** - Initializes with M equidistant levels

#### Source Coding
- **`Digital/Encoder/baseline_huffman_V2.m`** - Huffman encoder/decoder with Huffman tree construction
- **`Digital/Encoder/block_source_coding.m`** - Block-extension source coder (K-symbol blocks)

### Metrics & Utilities
- **`Metrics/MSE.m`** - Mean squared error
- **`Metrics/errorEnergy.m`** - Error energy (L2 norm squared)
- **`Math/integrate.m`** - Numerical integration via trapezoidal rule

### Phase 1: Testing & Analysis

- **`Testing Suites/sampling_analysis.m`** - Nyquist sampling, aliasing, reconstruction quality
- **`Testing Suites/fourier_analysis.m`** - Fourier series convergence analysis
- **`Testing Suites/quantization_analysis.m`** - Comparison of 2-level, Uniform, Lloyd-Max quantizers
- **`Testing Suites/coding_analysis.m`** - Huffman compression performance
- **`Testing Suites/coding_decoding_analysis.m`** - Full quantization + Huffman chain
- **`Testing Suites/block_coding_analysis.m`** - Block source coding with varying block sizes K

### Phase 2: Testing & Analysis

#### Modulation Demonstration
- **`Testing Suites/modulation_demo.m`** - BPSK/QPSK examples with constellation plots and noise analysis

#### BPSK Channel Simulation
- **`Testing Suites/simulate_bpsk_awgn.m`** - BPSK over AWGN channel, BER vs Eb/N0 plot
- **`Testing Suites/simulate_bpsk_coded_awgn.m`** - BPSK with repetition-3 coding over AWGN

#### QPSK Channel Simulation
- **`Testing Suites/simulate_qpsk_awgn.m`** - QPSK over AWGN channel, BER vs Eb/N0 plot
- **`Testing Suites/simulate_qpsk_coded_awgn.m`** - QPSK with repetition-3 coding over AWGN

#### ISI Channel & MLSE Equalization
- **`Testing Suites/simulate_qpsk_isi_viterbi.m`** - QPSK over ISI channel (y_i = 2x_i + x_{i-1} + n_i) with Viterbi MLSE
- **`Testing Suites/BER_estimation_qpsk_isi_mlse.m`** - Comprehensive QPSK ISI+MLSE BER analysis vs AWGN baseline

#### Verification
- **`Testing Suites/test_demod.m`** - Constellation and demodulation verification

### Deprecated
- **`util/decode_stream.m`** - Decodes bit stream using prefix-free code trie
- **`util/example1.m`** - Early testing example (reference)

---

## Execution

### Run Full Project (Recommended)
```matlab
>> main
```
Executes selected modules (Phase 1, Phase 2, or both based on uncommented lines).

### Run Individual Analysis
Make sure `constants.m` is added to `PATH`
```matlab
% Phase 1 Examples
>> sampling_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)
>> quantization_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)
>> fourier_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY)
>> coding_decoding_analysis(SIGNAL_DURATION, MESSAGE_FREQUENCY, CARRIER_FREQUENCY, TARGET_MSE, quantizer_levels)

% Phase 2 Examples
>> modulation_demo()
>> simulate_bpsk_awgn()
>> simulate_qpsk_awgn()
>> simulate_qpsk_isi_viterbi()
>> BER_estimation_qpsk_isi_mlse()
```

---

## Methods Implemented

### Phase 1: Signal Processing Chain
- **Sampling:** Nyquist-rate sampling with reconstruction via sinc interpolation
- **Fourier Analysis:** Fourier series with variable coefficients and period analysis
- **Quantization:**
  - Two-level quantizer (binary)
  - Uniform quantizer (M equally-spaced levels)
  - Lloyd-Max quantizer (optimal, iterative refinement)
- **Source Coding:**
  - Huffman compression (variable-length codes)
  - Block-extension coding (K-symbol dependencies)

### Phase 2: Digital Communications
- **Modulation:** BPSK (2-point), QPSK (4-point counter-clockwise)
- **Channels:**
  - AWGN (Additive White Gaussian Noise)
  - ISI (Intersymbol Interference): y_i = 2x_i + x_{i-1} + n_i
- **Equalization:** Viterbi MLSE (4-state trellis for QPSK with memory-1 ISI)
- **Channel Coding:** Repetition-3 (3× bit repetition)

### Metrics
- **Distortion:** MSE (mean squared error)
- **Throughput:** Bits per symbol
- **Latency:** Symbol-based and block-based delays
- **Performance:** BER (bit error rate) vs Eb/N0 (energy per bit to noise PSD ratio)
- **Error Energy:** L2 norm squared of approximation error

---

## Reproducibility

**No external toolboxes required.** All functions use core MATLAB functionality.

**MATLAB Version:** R2018b or later recommended (uses modern string handling, anonymous functions)

---

## Notes on AI Usage

We find LLMs create far better debug messages and comments than we can, and so LLMs have been used to provide comments for better readability. This documentation is also AI-generated, but checked by the project authors.