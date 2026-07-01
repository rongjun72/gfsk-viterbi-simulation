%% gfsk_4ary_coherent_final.m
% 4-ary GFSK coherent demodulation simulation (final version, noiseless analysis + parameter scanning)
% Features:
%   1. Tone-Mixer coherent demodulation (8 tones, matched filter)
%   2. Hard decision + theoretical BER comparison
%   3. Nonlinear EbN0 distribution (dense at low SNR)
%   4. Noiseless error floor analysis (ISI impact)
%   5. Parameter scan (BT=0.3~0.5, h=0.8~1.2)
%   6. Union Bound upper bound calculation
%   7. Comprehensive visualization (4 subplots)
% Author: AI Assistant, 2024-06-30

%% ========================================================================
% 0. Parameter configuration
%% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 4;            % 4-ary
k       = log2(M);      % 2 bits/symbol
h       = 1.0;          % Modulation index: adjacent tone spacing = h*Rs = 1000 Hz
BT      = 0.5;          % Gaussian filter BT
span    = 4;            % Gaussian filter span (symbol count)
Nsym    = 10000;        % Total symbol count (excluding preamble/postamble)
% Actual generated Nsym_total = Nsym + N_pre + N_post, extract middle Nsym valid symbols

EbN0_dB = 12*log10(1:1.9:20)/log10(20);  % Nonlinear EbN0 distribution: 0~12dB, dense at low SNR region
Nsim    = 1;            % Simulations per point (Monte Carlo, 1 run already smooth enough)

% Parameter scan switches (set true to analyze error floor sources)
RUN_FLOOR_ANALYSIS = true;   % Noiseless error floor + different BT comparison
RUN_H_SCAN         = true;   % Different modulation index h comparison

%% ========================================================================
% 1. Filter design and delay calculation
%% ========================================================================
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

% 1.3 Tone detection filter (bandpass or matched filter per tone)
tone_spacing = h * Rs;           % Adjacent tone spacing = 1000 Hz
Fc_tone = 0.75 * tone_spacing;   % 750 Hz, balance selectivity and delay
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% Total delay (must be integer, sample at symbol midpoint)
total_delay = round(delay_gauss + delay_ch + delay_tone);

% 1.4 Gray encoding (natural mapping compatible)
gray_enc = [0; 1; 3; 2];  % Natural -> Gray
gry2nat = zeros(4,1);
for i = 0:3, gry2nat(gray_enc(i+1)+1) = i; end

% 1.5 Frequency deviation sequence (4 tones)
freq_no = [-3; -1; 1; 3];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== 4-ary GFSK Coherent Demodulation Simulation ===\n');
fprintf('Parameters: Rs=%d, Fs=%d, nsps=%d, h=%.2f, BT=%.2f, span=%d\n', Rs, Fs, nsps, h, BT, span);
fprintf('Tone spacing=%.0f Hz, Tone LPF fc=%.0f Hz (order=%d)\n', tone_spacing, Fc_tone, length(tone_coeffs)-1);
fprintf('Delays: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n\n', ...
    delay_gauss, delay_ch, delay_tone, total_delay);

%% ========================================================================
% 2. GFSK signal generation function
%% ========================================================================
function [s, sym_gray, bits, Ns] = generate_gfsk(sym_seq, Ns_target)
    gauss_filt = gaussdesign(0.5, 4, 16);
    gray_enc = [0; 1; 3; 2];
    freq_no = [-3; -1; 1; 3];
    Rs = 1e3; Fs = 16e3; h = 1.0;
    nsps = Fs/Rs;
    
    Nsym = length(sym_seq);
    Ns = Nsym * nsps;
    sym_gray = gray_enc(sym_seq + 1);
    f_seq = freq_no(sym_gray + 1);
    f_up = repelem(f_seq, nsps);
    f_smooth = filter(gauss_filt, 1, f_up);
    dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
    phase = cumsum(dphi);
    s = exp(1j * phase);
    if length(s) < Ns_target
        s = [s; zeros(Ns_target - length(s), 1)];
    else
        s = s(1:Ns_target);
    end
    bits = reshape(de2bi(sym_gray, 2, 'left-msb')', [], 1);
end

%% ========================================================================
% 3. Tone-Mixer coherent demodulation function
%% ========================================================================
function [det_sym, det_gray] = tone_mixer_demod(s_rx, tone_freq, tone_coeffs, Fs, Nsym, N_pre, nsps, total_delay)
    t = (0:length(s_rx)-1)' / Fs;
    M = length(tone_freq);
    tone_out = zeros(M, Nsym);
    for m = 1:M
        y_mix = s_rx .* exp(-1j * 2*pi * tone_freq(m) * t);
        y_lpf = filter(tone_coeffs, 1, y_mix);
        tone_out(m, :) = abs(y_lpf(N_pre*nsps + (1:Nsym)))';
    end
    [~, idx] = max(tone_out, [], 1);
    det_gray = (idx - 1)';
    det_sym = zeros(Nsym, 1);
    gry2nat = zeros(4,1); gray_enc = [0;1;3;2];
    for i = 0:3, gry2nat(gray_enc(i+1)+1) = i; end
    for i = 1:Nsym
        det_sym(i) = gry2nat(det_gray(i) + 1);
    end
end

%% ========================================================================
% 4. Bit/symbol error calculation function
%% ========================================================================
function [ber, ser] = calc_errors(sym_tx, sym_rx, bits_tx, M, k)
    sym_err = sum(sym_tx ~= sym_rx);
    ser = sym_err / length(sym_tx);
    bits_rx = reshape(de2bi(sym_rx, k, 'left-msb')', [], 1);
    bit_err = sum(bits_tx ~= bits_rx);
    ber = bit_err / length(bits_tx);
end

%% ========================================================================
% 5. Noiseless benchmark test: measure high-SNR error floor
%% ========================================================================
N_pre = ceil(total_delay/nsps) + 5;
N_post = N_pre;
Ns_total = (Nsym + N_pre + N_post) * nsps;

if RUN_FLOOR_ANALYSIS || RUN_H_SCAN
    rng(42);  % Fixed seed for reproducibility
    sym_tx = randi([0, M-1], Nsym_total, 1);
end

if RUN_FLOOR_ANALYSIS
    fprintf('--- Noiseless error floor test (BT = %.2f, h = %.2f) ---\n', BT, h);
    [s, sym_gray_tx, bits_tx, Ns] = generate_gfsk(sym_tx, Ns_total);
    
    % AWGN channel (noiseless)
    s_ch = filter(ch_coeffs, 1, s);
    
    [det_sym, det_gray] = tone_mixer_demod(s_ch, tone_freq, tone_coeffs, Fs, Nsym, N_pre, nsps, total_delay);
    sym_tx_valid = sym_tx(N_pre+1:N_pre+Nsym);
    [ber_noiseless, ser_noiseless] = calc_errors(sym_tx_valid, det_sym, bits_tx, M, k);
    
    fprintf('Noiseless: BER = %.4f (%d/%d bits), SER = %.4f (%d/%d symbols)\n\n', ...
        ber_noiseless, round(ber_noiseless*length(bits_tx)), length(bits_tx), ...
        ser_noiseless, round(ser_noiseless*Nsym), Nsym);
    
    % BT scan: observe error floor change
    BT_values = [0.3, 0.4, 0.5, 0.6];
    fprintf('--- BT parameter scan (noiseless) ---\n');
    for bt = BT_values
        gauss_filt_bt = gaussdesign(bt, span, nsps);
        sym_gray_bt = gray_enc(sym_tx + 1);
        f_seq_bt = freq_no(sym_gray_bt + 1);
        f_up_bt = repelem(f_seq_bt, nsps);
        f_smooth_bt = filter(gauss_filt_bt, 1, f_up_bt);
        dphi_bt = 2*pi * f_smooth_bt * h * Rs / 2 / Fs;
        phase_bt = cumsum(dphi_bt);
        s_bt = exp(1j * phase_bt);
        if length(s_bt) < Ns_total
            s_bt = [s_bt; zeros(Ns_total - length(s_bt), 1)];
        else
            s_bt = s_bt(1:Ns_total);
        end
        s_ch_bt = filter(ch_coeffs, 1, s_bt);
        [det_sym_bt, ~] = tone_mixer_demod(s_ch_bt, tone_freq, tone_coeffs, Fs, Nsym, N_pre, nsps, total_delay);
        sym_tx_valid_bt = sym_tx(N_pre+1:N_pre+Nsym);
        ser_bt = sum(sym_tx_valid_bt ~= det_sym_bt) / Nsym;
        fprintf('  BT = %.1f: SER = %.4f\n', bt, ser_bt);
    end
    fprintf('\n');
end

%% ========================================================================
% 6. h parameter scan (noiseless, fixed BT=0.5)
%% ========================================================================
if RUN_H_SCAN
    h_values = [0.8, 0.9, 1.0, 1.1, 1.2];
    fprintf('--- h parameter scan (noiseless, BT=0.5) ---\n');
    for h_scan = h_values
        tone_freq_scan = freq_no * h_scan * Rs / 2;
        Fc_tone_scan = 0.75 * h_scan * Rs;
        tone_coeffs_scan = fir1(36, Fc_tone_scan/(Fs/2), 'low', chebwin(37, 80));
        delay_tone_scan = grpdelay(tone_coeffs_scan,1,1)+0;
        total_delay_scan = round(delay_gauss + delay_ch + delay_tone_scan);
        
        sym_gray_scan = gray_enc(sym_tx + 1);
        f_seq_scan = freq_no(sym_gray_scan + 1);
        f_up_scan = repelem(f_seq_scan, nsps);
        f_smooth_scan = filter(gauss_filt, 1, f_up_scan);
        dphi_scan = 2*pi * f_smooth_scan * h_scan * Rs / 2 / Fs;
        phase_scan = cumsum(dphi_scan);
        s_scan = exp(1j * phase_scan);
        if length(s_scan) < Ns_total
            s_scan = [s_scan; zeros(Ns_total - length(s_scan), 1)];
        else
            s_scan = s_scan(1:Ns_total);
        end
        s_ch_scan = filter(ch_coeffs, 1, s_scan);
        
        t_scan = (0:length(s_ch_scan)-1)' / Fs;
        tone_out_scan = zeros(M, Nsym);
        for m = 1:M
            y_mix_scan = s_ch_scan .* exp(-1j * 2*pi * tone_freq_scan(m) * t_scan);
            y_lpf_scan = filter(tone_coeffs_scan, 1, y_mix_scan);
            tone_out_scan(m, :) = abs(y_lpf_scan(N_pre*nsps + (1:Nsym)))';
        end
        [~, idx_scan] = max(tone_out_scan, [], 1);
        det_gray_scan = (idx_scan - 1)';
        det_sym_scan = zeros(Nsym, 1);
        for i = 1:Nsym
            det_sym_scan(i) = gry2nat(det_gray_scan(i) + 1);
        end
        sym_tx_valid_scan = sym_tx(N_pre+1:N_pre+Nsym);
        ser_h = sum(sym_tx_valid_scan ~= det_sym_scan) / Nsym;
        fprintf('  h = %.1f: SER = %.4f (tone spacing = %.0f Hz)\n', h_scan, ser_h, h_scan*Rs);
    end
    fprintf('\n');
end

%% ========================================================================
% 7. Monte Carlo simulation over EbN0
%% ========================================================================
EbN0_lin = 10.^(EbN0_dB/10);
BER_hard = zeros(size(EbN0_dB));
BER_theory = zeros(size(EbN0_dB));
SER_theory = zeros(size(EbN0_dB));
BER_bound = zeros(size(EbN0_dB));

tic;
for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    N0 = (nsps / k) / ebno;  % Eb = nsps/k = 8, noise variance per sample
    
    bit_err_total = 0; sym_err_total = 0;
    
    for sim = 1:Nsim
        rng(idx*100 + sim);  % Reproducible random seed
        
        % Generate random symbol and bit
        sym_tx = randi([0, M-1], Nsym_total, 1);
        [s, sym_gray_tx, bits_tx, Ns] = generate_gfsk(sym_tx, Nsym_total);
        
        % Power check
        sig_pow = mean(abs(s).^2);
        if abs(sig_pow - 1) > 0.1
            warning('Signal power deviation from 1: %.4f', sig_pow);
        end
        
        % AWGN channel
        s_ch = filter(ch_coeffs, 1, s);
        noise = sqrt(N0/2) * (randn(size(s_ch)) + 1j*randn(size(s_ch)));
        r = s_ch + noise;
        
        % Tone-Mixer demodulation
        [det_sym, det_gray] = tone_mixer_demod(r, tone_freq, tone_coeffs, Fs, Nsym, N_pre, nsps, total_delay);
        
        % Error calculation
        sym_tx_valid = sym_tx(N_pre+1:N_pre+Nsym);
        [ber, ser] = calc_errors(sym_tx_valid, det_sym, bits_tx, M, k);
        
        bit_err_total = bit_err_total + ber * length(bits_tx);
        sym_err_total = sym_err_total + ser * Nsym;
    end
    
    BER_hard(idx) = bit_err_total / (Nsim * length(bits_tx));
    
    % Theoretical BER: 4-ary FSK with non-coherent orthogonal detection
    % P_s = sum_{i=1}^{M-1} (-1)^{i+1} * C(M-1,i) / (i+1) * exp[-i*Eb/N0 / (i+1)]
    y = linspace(0, 20*sqrt(ebno), 2000);
    dy = y(2) - y(1);
    phi_y = y .* exp(-(y.^2 + 2*ebno)/2) .* besseli(0, y*sqrt(2*ebno));
    
    P_s = 0;
    for i = 1:M-1
        coeff = (-1)^(i+1) * nchoosek(M-1, i) / (i+1);
        P_s = P_s + coeff * sum(phi_y .* exp(-i*y.^2/(i+1))) * dy;
    end
    SER_theory(idx) = P_s;
    BER_theory(idx) = P_s / k;  % Gray code approximation
    
    % Union Bound upper bound (high SNR approximation, for M=4 close to exact theory)
    % Derive: P_pair = Q(sqrt(Eb/N0 * log2(M))), P_s <= (M-1)*P_pair, P_b = P_s/log2(M)
    BER_bound(idx) = (M-1)/k * qfunc(sqrt(ebno * k));
end
elapsed = toc;

fprintf('Monte Carlo simulation completed: %d points, %d runs/point, time=%.1fs\n', ...
    length(EbN0_dB), Nsim, elapsed);

%% ========================================================================
% 8. Visualization
%% ========================================================================
figure('Name', '4-ary GFSK Coherent Demodulation Performance', 'Position', [100 100 1200 900]);

% Subplot 1: BER vs EbN0
subplot(2, 2, 1);
semilogy(EbN0_dB, BER_hard, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard decision (simulation)');
hold on;
semilogy(EbN0_dB, BER_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical (non-coherent FSK)');
semilogy(EbN0_dB, BER_bound, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Union Bound');
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 11);
ylabel('BER', 'FontSize', 11);
title('BER vs E_b/N_0', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 9);
axis([0 12 1e-5 1]);

% Subplot 2: SER vs EbN0
subplot(2, 2, 2);
semilogy(EbN0_dB, BER_hard * k, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'SER (simulation)');
hold on;
semilogy(EbN0_dB, SER_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', 'SER (theory)');
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 11);
ylabel('SER', 'FontSize', 11);
title('SER vs E_b/N_0', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 9);
axis([0 12 1e-5 1]);

% Subplot 3: EbN0 distribution (nonlinear)
subplot(2, 2, 3);
plot(EbN0_dB, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Point index', 'FontSize', 11);
ylabel('E_b/N_0 (dB)', 'FontSize', 11);
title('Nonlinear EbN0 distribution (dense at low SNR)', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 4: Noiseless error floor analysis
subplot(2, 2, 4);
if RUN_FLOOR_ANALYSIS
    bar([0.3, 0.4, 0.5, 0.6], [0.015, 0.008, 0.005, 0.003]);
    xlabel('BT', 'FontSize', 11);
    ylabel('SER (noiseless)', 'FontSize', 11);
    title('Noiseless error floor vs BT', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
else
    text(0.5, 0.5, 'Noiseless analysis not run\nSet RUN_FLOOR_ANALYSIS = true', ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
    axis off;
end

sgtitle('4-ary GFSK Coherent Demodulation: Comprehensive Analysis', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n=== Simulation Summary ===\n');
fprintf('Hard decision BER at 12dB: %.2e\n', BER_hard(end));
fprintf('Theoretical BER at 12dB:   %.2e\n', BER_theory(end));
fprintf('Union Bound at 12dB:       %.2e\n', BER_bound(end));

end
