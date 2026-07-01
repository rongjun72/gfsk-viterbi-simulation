function gfsk_4ary_coherent_final()
% gfsk_4ary_coherent_final.m
% 4-ary GFSK coherent demodulation simulation — final archive version
% Includes: tone-mixer coherent detection, theoretical BER comparison, high-SNR error floor analysis
%
% Core architecture:
%   1. Tx: continuous-phase GFSK, gaussdesign(BT,span,sps) Gaussian pulse shaping
%   2. Channel: AWGN + 80dB out-of-band rejection channel filter
%   3. Rx: 4-branch tone-mixer coherent detection + Chebyshev window LPF + maximum magnitude decision
%   4. Analysis: theoretical M-ary orthogonal FSK coherent detection BER + noiseless error floor + BT parameter scan
%
% Key corrections (historical record):
%   - total_delay must include transmitter gauss_filt group delay (~16 samples),
%     otherwise sampling point systematically shifts by 1 symbol, causing BER≈0.5.
%   - Use dynamic N_pre/N_post = ceil(total_delay/nsps)+5 to ensure filter full steady-state.
%   - tone LPF fc = 0.75 * adjacent tone spacing, balancing adjacent channel isolation and ISI control.

%% ========================================================================
% 0. Configurable parameters
% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 4;            % 4-ary
k       = log2(M);      % 2 bits/symbol
h       = 1.0;          % Modulation index: adjacent tone spacing = h*Rs = 1000 Hz
BT      = 0.5;          % Gaussian filter BT
span    = 4;            % Gaussian filter span（Symbol count）
Nsym    = 10000;        % Total symbol count (excluding preamble/postamble)
% Actual generated Nsym_total = Nsym + N_pre + N_post, extract middle Nsym valid symbols

EbN0_dB = 12*log10(1:1.9:20)/log10(20);  % Nonlinear EbN0 distribution: 0~12dB, dense at low SNR region
Nsim    = 1;            % Simulations per point (Monte Carlo, 1 run already smooth enough)

% Parallel switch: requires Parallel Computing Toolbox
USE_PARFOR = true;      % true = parfor parallel acceleration for different EbN0 points; false = fallback to regular for

% Parameter scan switches (set true to analyze error floor sources)
RUN_FLOOR_ANALYSIS = true;   % Noiseless error floor + different BT comparison
RUN_H_SCAN         = true;   % Different modulation index h comparison

%% ========================================================================
% 1. Filter design and delay calculation
% ========================================================================
% 1.1 Gaussian frequency pulse (transmitter)
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

% 1.2 Channel filter: 80dB out-of-band rejection, lowpass FIR
% Optimization: Fp slightly larger than 4GFSK effective spectrum (±1750 Hz), reduce noise passage
Fp = 2.0e3;   Fs_stop = 2.8e3;  % Passband/stopband (Hz)
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
delay_ch = grpdelay(ch_filter.Coefficients,1,1)+0;
ch_coeffs = ch_filter.Coefficients;  % Extract coefficients, safe for parfor

% 1.3 Tone mixer lowpass filter: Chebyshev window, 36-tap, fc=0.75*tone spacing
tone_spacing = h * Rs;               % Adjacent tone spacing (Hz)
Fc_tone = 0.75 * tone_spacing;       % LPFCutoff frequency (Hz)
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% Total delay = tx Gaussian + channel + rx tone LPF
total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = ceil(total_delay/nsps) + 5;
Nsym_total = Nsym + N_pre + N_post;
Ns_total   = Nsym_total * nsps;      % Totalsample count

% Sampling instant: symbol midpoint + total delay compensation
sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay;
% Safety check
if sample_idx(1) < 1 || sample_idx(end) > Ns_total
    error('Sampling index out of bounds: total_delay=%d, N_pre=%d, firstIndex=%d, lastIndex=%d', ...
        total_delay, N_pre, sample_idx(1), sample_idx(end));
end

fprintf('=== 4-ary GFSK Coherent Demodulation ===\n');
fprintf('Parameters: Rs=%d, Fs=%d, nsps=%d, h=%.2f, BT=%.2f, span=%d\n', ...
    Rs, Fs, nsps, h, BT, span);
fprintf('Tone spacing=%.0f Hz, Tone LPF fc=%.0f Hz\n', tone_spacing, Fc_tone);
fprintf('Delays: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n', ...
    delay_gauss, delay_ch, delay_tone, total_delay);
fprintf('N_pre=%d, N_post=%d, Nsym_valid=%d, First sample=%d, Last sample=%d\n\n', ...
    N_pre, N_post, Nsym, sample_idx(1), sample_idx(end));

%% ========================================================================
% 2. Gray encoding/decoding mapping
% ========================================================================
% Natural binary -> Gray encoding
gray_enc = [0; 1; 3; 2];   % 00->0, 01->1, 11->3, 10->2
% Gray encoding -> natural binary
gry2nat = zeros(4,1);
for i = 0:3
    g = gray_enc(i+1);
    if g <= 3 && g >= 0 && g == round(g)
        gry2nat(g+1) = i;
    end
end
% Self-check
for i = 0:3
    if gry2nat(gray_enc(i+1)+1) ~= i
        error('Gray 映射不一致！');
    end
end

% Frequency indices (normalized, actual frequency = freq_no * h*Rs/2)
freq_no = [-3; -1; 1; 3];  % 4toneIndex

% 4 tones' actual frequencies (Hz), for mixing
tone_freq = freq_no * h * Rs / 2;

%% ========================================================================
% 3. Helper function: generate GFSK signal (noiseless)
% ========================================================================
    function [s, sym_gray, bits, Ns] = generate_gfsk(sym_seq, Ns_total)
        % sym_seq: 0..3 symbol sequence（LengthAny）
        Nsym_in = length(sym_seq);
        Ns_in = Nsym_in * nsps;
        
        % GrayEncode
        sym_gray_in = gray_enc(sym_seq + 1);
        % Frequency indexSequence
        f_seq = freq_no(sym_gray_in + 1);
        % UpsamplingAs pulse sequence
        f_up = repelem(f_seq, nsps);
        % Gaussian filter（Causal，Introducedelay_gauss）
        f_smooth = filter(gauss_filt, 1, f_up);
        % Phase integration（Continuous-phase）
        dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
        phase = cumsum(dphi);
        % Complex envelope signal
        s = exp(1j * phase);
        
        % EnsureLengthConsistent（May due toFilteringFilterGroup delaySlightly different）
        if length(s) < Ns_in
            s = [s; zeros(Ns_in - length(s), 1)];
        else
            s = s(1:Ns_in);
        end
        
        sym_gray = sym_gray_in;
        bits = zeros(Nsym_in*k, 1);
        for i = 1:Nsym_in
            nat = sym_seq(i);
            bits(2*i-1) = bitget(nat, 2);  % MSB
            bits(2*i)   = bitget(nat, 1);  % LSB
        end
        Ns = Ns_in;
    end

%% ========================================================================
% 4. Helper function: coherent detection (tone-mixer + LPF + sampling + decision)
% ========================================================================
    function [det_sym, det_gray, branch_metric] = detect_coherent(r, Nsym_in, sample_idx_in)
        % r: Received signal（AlreadyFilteringOrOriginal）
        % Nsym_in: Valid symbol count
        % sample_idx_in: Sampling index
        Ns_r = length(r);
        t = (0:Ns_r-1)' / Fs;
        
        branch_metric = zeros(M, Nsym_in);
        for m = 1:M
            % MixingTo basebandband
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            % LowpassFiltering（IsolateOther3tone）
            y_lpf = filter(tone_coeffs, 1, y_mix);
            % SymbolMidpointSample
            y_sample = y_lpf(sample_idx_in);
            % Branch metric = Magnitude
            branch_metric(m, :) = abs(y_sample).';  % 1 x Nsym
        end
        
        % Maximum magnitude decision -> GraySymbol（Force columnVectorAvoidDimensionNotMatch）
        [~, det_gray] = max(branch_metric, [], 1);
        det_gray = (det_gray(:) - 1);  % 0..3, ColumnVector Nsym x 1
        
        % Gray -> Natural binary
        det_sym = gry2nat(det_gray + 1);  % ColumnVector Nsym x 1
    end

if RUN_FLOOR_ANALYSIS || RUN_H_SCAN
    rng(42);  % Fixed seed for reproducibility
    sym_tx = randi([0, M-1], Nsym_total, 1);
end

%% ========================================================================
% 5. Noiseless benchmark test: measure high-SNR error floor
%% ========================================================================
if RUN_FLOOR_ANALYSIS
    fprintf('--- Noiseless error floor test (BT = %.2f, h = %.2f)---\n', BT, h);
    [s, ~, bits_tx, Ns] = generate_gfsk(sym_tx, Ns_total);
    
    % ChannelFiltering（Noiseless）
    r_ch = filter(ch_coeffs, 1, s);
    
    % Detect（UsingValidSymbolRange）
    [det_sym, det_gray, branch_metric] = detect_coherent(r_ch, Nsym, sample_idx);
    
    % CalculateError（OnlyCompareValidSymbol）
    sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
    sym_err = sum(det_sym ~= sym_tx_valid);
    
    % Re-CalculateBitError
    bits_tx_valid = zeros(Nsym*k, 1);
    for i = 1:Nsym
        nat = sym_tx_valid(i);
        bits_tx_valid(2*i-1) = bitget(nat, 2);
        bits_tx_valid(2*i)   = bitget(nat, 1);
    end
    bits_det_valid = zeros(Nsym*k, 1);
    for i = 1:Nsym
        nat = det_sym(i);
        bits_det_valid(2*i-1) = bitget(nat, 2);
        bits_det_valid(2*i)   = bitget(nat, 1);
    end
    bit_err = sum(bits_tx_valid ~= bits_det_valid);
    
    SER_floor = sym_err / Nsym;
    BER_floor = bit_err / (Nsym*k);
    fprintf('  SER_floor = %.4e, BER_floor = %.4e\n\n', SER_floor, BER_floor);
    
    % ---- Different BT error floorScan ----
    BT_scan = [0.3, 0.5, 1.0];
    SER_floor_BT = zeros(size(BT_scan));
    BER_floor_BT = zeros(size(BT_scan));
    fprintf('--- 不同 BT 的误差地板扫描（h = %.2f）---\n', h);
    for bt_idx = 1:length(BT_scan)
        BT_tmp = BT_scan(bt_idx);
        gauss_tmp = gaussdesign(BT_tmp, span, nsps);
        
        % Regenerate signal（OnlyBTChange）
        sym_gray_tmp = gray_enc(sym_tx + 1);
        f_seq_tmp = freq_no(sym_gray_tmp + 1);
        f_up_tmp = repelem(f_seq_tmp, nsps);
        f_smooth_tmp = filter(gauss_tmp, 1, f_up_tmp);
        dphi_tmp = 2*pi * f_smooth_tmp * h * Rs / 2 / Fs;
        phase_tmp = cumsum(dphi_tmp);
        s_tmp = exp(1j * phase_tmp);
        s_tmp = s_tmp(1:Ns_total);
        
        r_ch_tmp = filter(ch_coeffs, 1, s_tmp);
        [det_sym_tmp, ~, ~] = detect_coherent(r_ch_tmp, Nsym, sample_idx);
        sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
        SER_floor_BT(bt_idx) = sum(det_sym_tmp ~= sym_tx_valid) / Nsym;
        % BitError
        bits_err_bt = 0;
        for i = 1:Nsym
            bits_err_bt = bits_err_bt + (bitget(sym_tx_valid(i),2) ~= bitget(det_sym_tmp(i),2)) ...
                + (bitget(sym_tx_valid(i),1) ~= bitget(det_sym_tmp(i),1));
        end
        BER_floor_BT(bt_idx) = bits_err_bt / (Nsym*k);
        
        fprintf('  BT=%.1f: SER_floor=%.4e, BER_floor=%.4e\n', ...
            BT_tmp, SER_floor_BT(bt_idx), BER_floor_BT(bt_idx));
    end
    fprintf('\n');
end

%% ========================================================================
% 6. Error floor scan for different modulation index h
%% ========================================================================
if RUN_H_SCAN
    h_scan = [0.5, 1.0];
    SER_floor_h = zeros(size(h_scan));
    BER_floor_h = zeros(size(h_scan));
    fprintf('--- 不同 h 的误差地板扫描（BT = %.2f）---\n', BT);
    for h_idx = 1:length(h_scan)
        h_tmp = h_scan(h_idx);
        tone_spacing_tmp = h_tmp * Rs;
        tone_freq_tmp = freq_no * h_tmp * Rs / 2;
        Fc_tone_tmp = 0.75 * tone_spacing_tmp;
        tone_coeffs_tmp = fir1(36, Fc_tone_tmp/(Fs/2), 'low', chebwin(37, 80));
        
        % Delay re-Calculate（toneFilteringFilter unchanged，But frequency changed，LengthUnchanged，delayUnchanged）
        % Onlytone_freqChanged，tone_coeffsRedesign
        
        % Generate signal（h_tmp）
        sym_gray_tmp = gray_enc(sym_tx + 1);
        f_sel_tmp = freq_no(sym_gray_tmp + 1);  % SelectofNormalizationFrequency index
        f_up_tmp = repelem(f_sel_tmp, nsps);
        f_smooth_tmp = filter(gauss_filt, 1, f_up_tmp);
        dphi_tmp = 2*pi * f_smooth_tmp * h_tmp * Rs / 2 / Fs;
        phase_tmp = cumsum(dphi_tmp);
        s_tmp = exp(1j * phase_tmp);
        s_tmp = s_tmp(1:Ns_total);
        
        r_ch_tmp = filter(ch_coeffs, 1, s_tmp);
        
        % Use newtone_freqDetect
        Ns_r = length(r_ch_tmp);
        t = (0:Ns_r-1)' / Fs;
        branch_m_h = zeros(M, Nsym);
        for m = 1:M
            y_mix = r_ch_tmp .* exp(-1j * 2*pi * tone_freq_tmp(m) * t);
            y_lpf = filter(tone_coeffs_tmp, 1, y_mix);
            y_sample = y_lpf(sample_idx);
            branch_m_h(m, :) = abs(y_sample).';
        end
        [~, det_gray_h] = max(branch_m_h, [], 1);
        det_gray_h = det_gray_h(:) - 1;
        det_sym_h = gry2nat(det_gray_h + 1);
        
        sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
        SER_floor_h(h_idx) = sum(det_sym_h ~= sym_tx_valid) / Nsym;
        bits_err_h = 0;
        for i = 1:Nsym
            bits_err_h = bits_err_h + (bitget(sym_tx_valid(i),2) ~= bitget(det_sym_h(i),2)) ...
                + (bitget(sym_tx_valid(i),1) ~= bitget(det_sym_h(i),1));
        end
        BER_floor_h(h_idx) = bits_err_h / (Nsym*k);
        
        fprintf('  h=%.1f: SER_floor=%.4e, BER_floor=%.4e\n', ...
            h_tmp, SER_floor_h(h_idx), BER_floor_h(h_idx));
    end
    fprintf('\n');
end

%% ========================================================================
% 7. Main simulation: Eb/N0 sweep
%% ========================================================================
EbN0_lin = 10.^(EbN0_dB/10);
BER_sim = zeros(size(EbN0_dB));
SER_sim = zeros(size(EbN0_dB));

if USE_PARFOR
    % CheckAnd start parallel pool（Need Parallel Computing Toolbox）
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local');
        pool = gcp('nocreate');
    end
    fprintf('Using parfor with %d workers on %d EbN0 points x %d sims each\n', ...
        pool.NumWorkers, length(EbN0_dB), Nsim);
    
    t_start = tic;
    parfor idx = 1:length(EbN0_dB)
        ebno = EbN0_lin(idx);
        N0 = 8 / ebno;
        noise_var = N0;
        
        bit_err_local = 0;
        sym_err_local = 0;
        
        for sim = 1:Nsim
            rng(idx*100 + sim, 'twister');
            sym_tx = randi([0, M-1], Nsym_total, 1);
            
            % === Inline generate_gfsk（parfor Does not support nestingFunction）===
            Nsym_in = length(sym_tx);
            Ns_in = Nsym_in * nsps;
            sym_gray_in = gray_enc(sym_tx + 1);
            f_seq = freq_no(sym_gray_in + 1);
            f_up = repelem(f_seq, nsps);
            f_smooth = filter(gauss_filt, 1, f_up);
            dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
            phase = cumsum(dphi);
            s = exp(1j * phase);
            if length(s) < Ns_in
                s = [s; zeros(Ns_in - length(s), 1)];
            else
                s = s(1:Ns_in);
            end
            
            sig_pow = mean(abs(s).^2);
            if abs(sig_pow - 1) > 0.1
                warning('Signal power deviation: %.3f', sig_pow);
            end
            
            noise = sqrt(noise_var/2) * (randn(Ns_total, 1) + 1j*randn(Ns_total, 1));
            r = s + noise;
            r_ch = filter(ch_coeffs, 1, r);
            
            % === Inline detect_coherent（parfor Does not support nestingFunction）===
            Ns_r = length(r_ch);
            t = (0:Ns_r-1)' / Fs;
            branch_metric = zeros(M, Nsym);
            for m = 1:M
                y_mix = r_ch .* exp(-1j * 2*pi * tone_freq(m) * t);
                y_lpf = filter(tone_coeffs, 1, y_mix);
                y_sample = y_lpf(sample_idx);
                branch_metric(m, :) = abs(y_sample).';
            end
            [~, det_gray] = max(branch_metric, [], 1);
            det_gray = det_gray(:) - 1;
            det_sym = gry2nat(det_gray + 1);
            
            sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
            sym_err_local = sym_err_local + sum(det_sym ~= sym_tx_valid);
            
            for i = 1:Nsym
                bit_err_local = bit_err_local + ...
                    (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym(i), 2)) + ...
                    (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym(i), 1));
            end
        end
        
        BER_sim(idx) = bit_err_local / (Nsym * k * Nsim);
        SER_sim(idx) = sym_err_local / (Nsym * Nsim);
    end
    t_total = toc(t_start);
    fprintf('Total parfor time: %.2f s\n', t_total);
else
    % Serial fallback（With gfsk_4ary_coherent_final.m Behavior consistent）
    tic;
    for idx = 1:length(EbN0_dB)
        ebno = EbN0_lin(idx);
        N0 = 8 / ebno;
        noise_var = N0;
        
        bit_err_total = 0;
        sym_err_total = 0;
        
        for sim = 1:Nsim
            rng(idx*100 + sim);
            sym_tx = randi([0, M-1], Nsym_total, 1);
            [s, ~, ~, ~] = generate_gfsk(sym_tx, Ns_total);
            
            sig_pow = mean(abs(s).^2);
            if abs(sig_pow - 1) > 0.1
                warning('Signal power deviation: %.3f', sig_pow);
            end
            
            noise = sqrt(noise_var/2) * (randn(Ns_total, 1) + 1j*randn(Ns_total, 1));
            r = s + noise;
            r_ch = filter(ch_coeffs, 1, r);
            
            [det_sym, ~, ~] = detect_coherent(r_ch, Nsym, sample_idx);
            
            sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
            sym_err_total = sym_err_total + sum(det_sym ~= sym_tx_valid);
            
            for i = 1:Nsym
                bit_err_total = bit_err_total + ...
                    (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym(i), 2)) + ...
                    (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym(i), 1));
            end
        end
        
        BER_sim(idx) = bit_err_total / (Nsym * k * Nsim);
        SER_sim(idx) = sym_err_total / (Nsym * Nsim);
        
        fprintf('Eb/N0=%5.1f dB | BER=%.4e | SER=%.4e | Time=%.2f s\n', ...
            EbN0_dB(idx), BER_sim(idx), SER_sim(idx), toc);
    end
end

%% ========================================================================
% 8. Theoretical BER calculation: M-ary orthogonal FSK coherent detection (Gray coding)
%% ========================================================================
% Exact formula：P_s = 1 - integral_{-inf}^{inf} phi(y) * [1-Q(y+sqrt(2*ebno*log2(M)))]^(M-1) dy
% Where phi(y) = normpdf(y), Q(y) = qfunc(y)
% P_b ≈ P_s / log2(M) （HighSNRApproximate）

BER_theory = zeros(size(EbN0_dB));
SER_theory = zeros(size(EbN0_dB));

for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    % UsingNumerical integrationCalculateExact SER
    % IntegrandFunction：phi(y) * [1 - Q(y + sqrt(2*ebno*log2(M)))]^(M-1)
    sqrt_term = sqrt(2 * ebno * log2(M));
    
    % IntegrateRangeTake [-5*sqrt_term-10, 5] Usually sufficient
    y_min = -5*sqrt_term - 10;
    y_max = 10;
    if y_min < -50, y_min = -50; end
    
    % Using Simpson Or trapezoidal rule（Avoid dependency integral Toolbox）
    Ny = 2000;
    y = linspace(y_min, y_max, Ny);
    dy = y(2) - y(1);
    phi_y = (1/sqrt(2*pi)) * exp(-y.^2/2);
    Q_term = qfunc(y + sqrt_term);
    integrand = phi_y .* (1 - Q_term).^(M-1);
    P_s = 1 - trapz(y, integrand);
    
    SER_theory(idx) = P_s;
    BER_theory(idx) = P_s / log2(M);  % Gray code approximation
end

% Union bound upper bound (high SNR approximation, for M=4 close to exact theory)
% Derive：P_pair = Q(sqrt(Eb/N0 * log2(M))), P_s <= (M-1)*P_pair, P_b = P_s/log2(M)
BER_bound = (M-1)/log2(M) * qfunc(sqrt(EbN0_lin * log2(M)));

%% ========================================================================
% 9. Visualization
%% ========================================================================

% Figure 1: BER curves (simulation vs theory)
figure('Name', 'BER/SER Performance', 'Position', [100 100 800 600]);
semilogy(EbN0_dB, BER_sim, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Simulated BER');
hold on;
semilogy(EbN0_dB, SER_sim, 'rs--', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Simulated SER');
semilogy(EbN0_dB, BER_theory, 'g-', 'LineWidth', 1, 'DisplayName', 'Theory BER (orthog. MFSK coherent)');
semilogy(EbN0_dB, BER_bound, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Union Bound');

% Mark noiseless floor
if RUN_FLOOR_ANALYSIS
    xline(max(EbN0_dB)+2, 'k--', 'DisplayName', 'Noiseless floor region');
    text(max(EbN0_dB)+1, BER_floor, sprintf('Floor: %.2e', BER_floor), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate / Symbol Error Rate');
legend('Location', 'southwest');
title(sprintf('4-ary GFSK Coherent Detection (h=%.1f, BT=%.1f, nsps=%d)', h, BT, nsps));
axis([min(EbN0_dB) max(EbN0_dB) 1e-4 1]);

% Figure 2: Error floor analysis (BT and h scan)
if RUN_FLOOR_ANALYSIS && RUN_H_SCAN
    figure('Name', 'Error Floor Analysis', 'Position', [150 150 900 400]);
    
    subplot(1,2,1);
    bar(BT_scan, BER_floor_BT);
    set(gca, 'YScale', 'log');
    xlabel('BT');
    ylabel('BER Floor');
    title(sprintf('Error Floor vs BT (h=%.1f)', h));
    grid on;
    for i = 1:length(BT_scan)
        text(i, BER_floor_BT(i)*1.5, sprintf('%.2e', BER_floor_BT(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    
    subplot(1,2,2);
    bar(h_scan, BER_floor_h);
    set(gca, 'YScale', 'log');
    xlabel('h');
    ylabel('BER Floor');
    title(sprintf('Error Floor vs h (BT=%.1f)', BT));
    grid on;
    for i = 1:length(h_scan)
        text(i, BER_floor_h(i)*1.5, sprintf('%.2e', BER_floor_h(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

% Figure 3: Single waveform illustration (last EbN0 point, if Nsim>=1)
figure('Name', 'Signal Waveform & Detection', 'Position', [200 200 1000 700]);

% Regenerate a noiseless/low-noise example signal for plotting
rng(123);
sym_demo = randi([0, M-1], Nsym_total, 1);
[s_demo, ~, ~, ~] = generate_gfsk(sym_demo, Ns_total);
r_demo = filter(ch_coeffs, 1, s_demo);
[det_sym_demo, det_gray_demo, branch_metric_demo] = detect_coherent(r_demo, Nsym, sample_idx);

% Time axis
t_demo = (0:length(r_demo)-1)/Fs;
t_sym = (N_pre + (0:Nsym-1)) / Rs;

% 3.1 Transmitted signal spectrum
subplot(3,2,1);
[pxx, f] = periodogram(s_demo, hamming(length(s_demo)), 4096, Fs, 'centered');
plot(f/1e3, 10*log10(pxx));
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('Transmitted 4-GFSK Spectrum');
grid on; xlim([-4 4]);

% 3.2 Channel filter response
subplot(3,2,2);
[Hch, fch] = freqz(ch_filter.Coefficients, 1, 4096, Fs);
plot(fch/1e3, 20*log10(abs(Hch)));
hold on;
for m = 1:M
    xline(tone_freq(m)/1e3, 'r--', sprintf('f_%d', m));
end
xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
title('Channel Filter Response & Tone Frequencies');
grid on; xlim([-4 4]); ylim([-100 5]);

% 3.3 4-branch output (time domain)
subplot(3,2,3);
plot(t_demo*1e3, real(r_demo), 'b', 'LineWidth', 0.5);
hold on;
plot(t_demo(sample_idx)*1e3, real(r_demo(sample_idx)), 'ro', 'MarkerSize', 6);
xlabel('Time (ms)'); ylabel('Real(r(t))');
title('Received Signal & Sampling Points');
grid on; xlim([0 5]);

% 3.4 4-branch metric values
subplot(3,2,4);
plot(t_sym*1e3, branch_metric_demo.', 'LineWidth', 1.5);
hold on;
plot(t_sym*1e3, det_gray_demo+1, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
xlabel('Time (ms)'); ylabel('Branch Metric');
legend('Branch 0', 'Branch 1', 'Branch 2', 'Branch 3', 'Decision', 'Location', 'best');
title('Tone-Mixer Branch Metrics at Symbol Midpoints');
grid on; xlim([0 5]);

% 3.5 Transmitted vs detected symbol comparison
subplot(3,2,5);
stem(1:50, sym_demo(N_pre+1:N_pre+50), 'b', 'LineWidth', 1.5, 'DisplayName', 'Tx');
hold on;
stem(1:50, det_sym_demo(1:50), 'r--', 'LineWidth', 1, 'DisplayName', 'Detected');
xlabel('Symbol Index'); ylabel('Symbol Value');
title('Transmitted vs Detected Symbols (First 50)');
legend('Location', 'best'); grid on;

% 3.6 Constellation/eye diagram illustration (I-Q trajectory)
subplot(3,2,6);
plot(real(r_demo), imag(r_demo), 'b.', 'MarkerSize', 3);
hold on;
plot(real(r_demo(sample_idx)), imag(r_demo(sample_idx)), 'ro', 'MarkerSize', 6);
xlabel('In-Phase'); ylabel('Quadrature');
title('Constellation Diagram at Sampling Points');
axis equal; grid on;

%% ========================================================================
% 10. Results summary output
%% ========================================================================
fprintf('\n========== RESULT SUMMARY ==========\n');
fprintf('Eb/N0(dB) |  Sim BER  |  Sim SER  | Theory BER | Union Bound\n');
fprintf('----------|-----------|-----------|------------|------------\n');
for idx = 1:length(EbN0_dB)
    fprintf('%7.1f   | %.4e | %.4e |  %.4e  | %.4e\n', ...
        EbN0_dB(idx), BER_sim(idx), SER_sim(idx), BER_theory(idx), BER_bound(idx));
end

if RUN_FLOOR_ANALYSIS
    fprintf('\n----- Error Floor Analysis (Noiseless) -----\n');
    fprintf('Default (BT=%.2f, h=%.1f): BER_floor = %.4e, SER_floor = %.4e\n', ...
        BT, h, BER_floor, SER_floor);
    fprintf('\nBT Scan:\n');
    for i = 1:length(BT_scan)
        fprintf('  BT=%.1f: BER_floor = %.4e, SER_floor = %.4e\n', ...
            BT_scan(i), BER_floor_BT(i), SER_floor_BT(i));
    end
end

if RUN_H_SCAN
    fprintf('\nh Scan:\n');
    for i = 1:length(h_scan)
        fprintf('  h=%.1f: BER_floor = %.4e, SER_floor = %.4e\n', ...
            h_scan(i), BER_floor_h(i), SER_floor_h(i));
    end
end

fprintf('\nSimulation complete.\n');

end
