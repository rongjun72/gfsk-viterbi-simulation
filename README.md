# GFSK Viterbi Simulation

MATLAB simulation framework for coherent demodulation of 4-ary and 8-ary GFSK (Gaussian Frequency Shift Keying) with ISI-aware Viterbi sequence detection.

## Project Overview

This project implements a complete end-to-end MATLAB simulation system for GFSK coherent demodulation, from signal generation through AWGN channel transmission to receiver detection. The core innovation is an **ISI-aware Viterbi sequence detector** designed to combat the inter-symbol interference introduced by Gaussian pulse shaping.

### Key Features

- **4-ary and 8-ary GFSK modulation** with continuous-phase generation
- **Coherent tone-mixer detection** with Chebyshev-windowed LPF
- **ISI-aware Viterbi decoding** with pre-computed reference templates
- **Extensive parameter optimization** tools (h, LPF order, cutoff frequency)
- **Sample phase analysis** and per-branch optimal phase investigation
- **Python MLSE tutorial** with BPSK + ISI channel examples

## System Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Symbol Rate (Rs) | 1 kHz | Symbol duration 1 ms |
| Sampling Rate (Fs) | 16 kHz | 16× oversampling |
| Modulation Index (h) | 1.0 (default) | Adjustable: 0.9 - 1.2 |
| Gaussian BT | 0.5 | Fixed constraint |
| Gaussian Span | 4 | Symbols |
| Tone LPF Order | 24 (default) | Even orders recommended |
| Tone LPF fc | 0.75 × tone_spacing | Adjustable |

### 4-ary Configuration
- 4 tones: [-1500, -500, 500, 1500] Hz
- Gray code: [0, 1, 3, 2]
- Channel filter: Fp=2.0 kHz, Fs_stop=2.8 kHz

### 8-ary Configuration
- 8 tones: [-3500, -2500, -1500, -500, 500, 1500, 2500, 3500] Hz
- Gray code: [0, 1, 3, 2, 6, 7, 5, 4]
- Channel filter: Fp=4.5-5.0 kHz, Fs_stop=5.5-6.0 kHz

## File Structure

### Main Simulation Scripts (MATLAB)

| File | Description | Status |
|------|-------------|--------|
| `gfsk_4ary_coherent_final.m` | 4-ary hard-decision baseline + theoretical BER | Final |
| `gfsk_4ary_viterbi_isi.m` | **4-state ISI-aware Viterbi** (core result) | Final |
| `gfsk_4ary_viterbi_64state.m` | 64-state Viterbi (3-symbol memory, experimental) | Experimental |
| `gfsk_4ary_coherent_parfor.m` | Parallel hard-decision version | Final |
| `gfsk_4ary_viterbi.m` | Basic Viterbi (no ISI awareness, historical) | Historical |
| `gfsk_8ary_coherent_final.m` | 8-ary hard-decision baseline | Final |
| `gfsk_8ary_viterbi_isi.m` | **8-state ISI-aware Viterbi** (core result) | Final |
| `gfsk_8ary_viterbi_64state.m` | 64-state 8-ary Viterbi (experimental) | Experimental |

### Analysis & Optimization Scripts

| File | Description |
|------|-------------|
| `plot_viterbi_trellis.m` | Viterbi trellis diagram visualization (3 figures) |
| `analyze_symbol_internal_metrics.m` | 4-ary symbol-internal branch metric analysis |
| `test_tone_lpf_order.m` | Tone LPF order sweep (10-50) with BER comparison |
| `analyze_8gfsk_tone_metrics_intrasync.m` | 8-ary intra-symbol tone metric distribution |
| `analyze_per_branch_optimal_phase.m` | Per-branch optimal phase test (8-ary) |
| `analyze_per_branch_optimal_phase_4ary.m` | Per-branch optimal phase test (4-ary) |
| `scan_sample_phase_8gfsk.m` | Sample phase offset scan (-3 to +3) |
| `optimize_8gfsk_h_lpf.m` | Joint optimization of h, LPF order, and fc factor |

### Documentation

| File | Description |
|------|-------------|
| `gfsk_4ary_simulation_report.md` | Complete simulation report with 6 figures |
| `MLSE_Tutorial_temp.md` | Maximum Likelihood Sequence Estimation tutorial |
| `Trellis_temp.md` | Trellis state transition detailed explanation |

### Python MLSE Tutorial (`MLSE_Viterbi_Tutorial/`)

| File | Description |
|------|-------------|
| `MLSE_Viterbi_Tutorial/MLSE教程.md` | MLSE tutorial (Chinese) |
| `MLSE_Viterbi_Tutorial/Trellis状态转移详解.md` | Trellis explanation (Chinese) |
| `MLSE_Viterbi_Tutorial/mlse_viterbi_demo.py` | BPSK MLSE vs MMSE BER comparison |
| `MLSE_Viterbi_Tutorial/sliding_window_viterbi.py` | Real-time sliding-window Viterbi |
| `MLSE_Viterbi_Tutorial/viterbi_forward_trace.py` | Step-by-step Viterbi forward trace |

## Key Results

### 4-ary GFSK
- Hard-decision: BER ≈ 0% noiseless (perfect tone spacing)
- Viterbi gain: ~3-5 dB at high SNR, reduces error floor by ~1 order of magnitude

### 8-ary GFSK (h=1.0, order=24, fc=0.75)
- Hard-decision: SER ≈ 16.25% noiseless (tone crowding + Gaussian ISI)
- Viterbi: SER ≈ 1.08% noiseless (significant ISI mitigation)

### 8-ary Optimization (h=1.2, order=20, fc=1.0)
- Best SER: **3.5%** noiseless (vs 6.1% at h=1.0, order=20, fc=1.0)
- Key insight: **lower delay dominates** over selectivity; wide fc + low order = better sampling position

### Critical Design Notes

1. **Total delay must include all three filters**: `delay_gauss + delay_ch + delay_tone`
   - Missing `delay_gauss` causes BER ≈ 0.5 (systematic failure)

2. **Tone LPF order must be even**
   - Odd orders have half-integer `grpdelay` → `round(total_delay)` creates 0.5-sample offset
   - Even orders 10-44 pass noiseless check reliably

3. **Viterbi reference templates must include channel filter**
   - Without `ch_filter`, templates mismatch real signal by amplitude/frequency response

4. **Preamble must be fixed `0`**
   - Random preamble causes Viterbi 50% error (reference assumes prev=0)

5. **Sample phase matters**
   - Negative delta (-1 to -3): Viterbi works (ISI from past only, predictable)
   - Positive delta (+1 to +3): Viterbi fails (future leakage, unpredictable)

## Requirements

- MATLAB with Signal Processing Toolbox (`gaussdesign`, `designfilt`, `fir1`, `chebwin`, `grpdelay`)
- Python 3 with numpy and matplotlib (for tutorial scripts)
- Optional: Parallel Computing Toolbox (`parfor` in `gfsk_4ary_coherent_parfor.m`)

## Usage

```matlab
% Run 4-ary simulation with ISI-aware Viterbi
gfsk_4ary_viterbi_isi

% Run 8-ary simulation
gfsk_8ary_viterbi_isi

% Optimize 8-ary parameters
optimize_8gfsk_h_lpf

% Visualize Viterbi trellis
plot_viterbi_trellis
```

## License

This project is created for research and educational purposes.

## Acknowledgments

- G. D. Forney Jr., "The Viterbi Algorithm" (1973)
- Proakis & Salehi, *Digital Communications*
- MATLAB Signal Processing Toolbox documentation
