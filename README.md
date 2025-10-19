# Digital Communication Chain (Phase 1)

**Authors:** Khalil El Kaaki, Mouhammad Kandakji
**Course:** EECE 442: Communication Systems
**Tags:** Sampling, Fourier Series, Quantization, Lossless Source Coding Chain
**Institution:** American University of Beirut

---

## Overview

This project implements a basic signal processing chain that includes analog signal sampling, quantization, lossless compression, and reconstruction.
The main goal is to study trade-offs between quantization rate, distortion (MSE), and source statistics.

---

## Repository Structure

```
├── main.m       # Full-chain experiment (sampling → quantization → Huffman)
├── util/
|   ├── baseline_huffman_V2.m  # Huffman encoder/decoder
|   ├── block_source_coding.m  # Implements a block source coder of size K.
|   ├── decodeStream.m         # Decodes a bitstream to its quantization levels
|   ├── errEn.m                # Calculates the error energy
|   ├── example1.m             # Was used for early testing (deprecated)
|   ├── exampleSpeechWave.m    # Generates an AM-modulated wave
|   ├── fcoef.m                # Helper for fs.m
|   ├── fsConstT.m             # Plots fourier series for a fixed period T
|   ├── fsConstN.m             # Plots fourier series for a fixed # coefficients N
|   ├── fs.m                   # Calculates the fourier series of a signal
|   ├── lloydMaxInit.m         # Helper for lloydMax.m
|   ├── lloydMax.m             # Implements a Lloyd-Max Quantizer
|   ├── MSE.m                  # Calculates a sample's MSE
|   ├── partition.m            # Helper for lloydMax.m
|   ├── quan.m                 # Helper for UniformQuan.m
|   ├── reconstruct.m          # Reconstruct signal from samples
|   ├── sample.m               # Analog Sampler
|   ├── sinc.m                 # sinc(t/T) function definition
|   ├── twoLvlQuan.m           # Two-level Quantizer
|   ├── UniformQuan.m          # Uniform Quantizer
│   └── plotting scripts       # Visualizations
├── sampling_analysis.m    # Samples a signal and reconstructs it at different frequencies
├── fourierAnalysis.m      # Generate approximate signals from the Fourier Series
├── quantization_analysis.m# Tests all the quantizers and plots figures
├── block_coding_analysis.m# Block coding and performance analysis (deprecated)
├── coding_analysis.m      # Baseline coding and performance analysis (deprecated)
└── README.md
```

---

## Execution

From the MATLAB console, run `main.m`:

```matlab
>> main
```

It generates:

1. Plots of all the signals, their samples, approximates, reconstruction, and the variation of MSE, throughput with the variation in parameters.
2. Summary tables printed to console.

---

## Methods Implemented

* **Signal Sampling:** Nyquist-rate sampling of an AM-like waveform.
* **Quantization:**
  * Two-level quantizer.
  * Uniform quantizer.
  * Lloyd–Max quantizer with target MSE.
* **Lossless Coding:**
  * Huffman compression.
  * Block-coding extension.
* **Metrics:**
  * Total bits
  * Throughput (bits/symbol)
  * Latency (symbols)
  * End-to-end distortion (MSE)

---

## Reproducibility

No external toolboxes required.

---
