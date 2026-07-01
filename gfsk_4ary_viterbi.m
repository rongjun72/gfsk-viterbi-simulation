function gfsk_4ary_viterbi()
% gfsk_4ary_viterbi.m
% 4-ary GFSK coherent demodulation + Viterbi sequence detection simulation
%
% Architecture:
%   1. Tx: continuous-phase GFSK, gaussdesign(BT,span,sps) Gaussian pulse shaping
%   2. Channel: AWGN + 80dB out-of-band rejection channel filter
%   3. Receiver frontend：4-branch tone-mixer coherent detection + Chebyshev window LPF
%   4. Rx back-end: 4-state Viterbi soft-decision sequence detection (against ISI)
%   5. Comparison: per-symbol hard decision vs Viterbi sequence detection
%
% Key design:
%   - Viterbi state = current symbol candidate frequency (0..3)
%   - Branch metric = 4 magnitude values from tone-mixer at each symbol midpoint
%   - Path normalization prevents numerical overflow
%   - Traceback length = full frame (Nsym), standard survivor path traceback

%% ========================================================================
% 0. Configurable parameters
%% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 4;            % 4-ary
k       = log2(M);      % 2 bits/symbol
h       = 1.0;          % Modulation index: adjacent tone spacing = h*Rs = 1000 Hz
BT      = 0.5;          % Gaussian filter BT
span    = 4;            % Gaussian filter span（Symbol count）
Nsym    = 10000;        % ValidSymbol count

EbN0_dB = 12*log10(1:1.9:20)/log10(20);  % Nonlinear EbN0 distribution: 0~12dB
Nsim    = 1;            % Simulations per point

%% ========================================================================
% 1. Filter design and delay calculation
%% ========================================================================
% 1.1 Gaussian frequency pulse (transmitter)
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

% 1.2 Channel filter: 80dB out-of-band rejection, lowpass FIR
% Optimization: Fp slightly larger than 4GFSK effective spectrum (±1750 Hz), reduce noise passage
Fp = 2.0e3;   Fs_stop = 2.8e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
delay_ch = grpdelay(ch_filter.Coefficients,1,1)+0;

% 1.3 Tone mixer lowpass filter: Chebyshev window, 36-tap, fc=0.75*tone spacing
tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% Total delay
total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = ceil(total_delay/nsps) + 5;
Nsym_total = Nsym + N_pre + N_post;
Ns_total   = Nsym_total * nsps;

% Sampling instant: symbol midpoint + total delay compensation
sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay;
if sample_idx(1) < 1 || sample_idx(end) > Ns_total
    error('Sampling index out of bounds: total_delay=%d, N_pre=%d, firstIndex=%d, lastIndex=%d', ...
        total_delay, N_pre, sample_idx(1), sample_idx(end));
end

fprintf('=== 4-ary GFSK Coherent + Viterbi Sequence Detection ===\n');
fprintf('Parameters: Rs=%d, Fs=%d, nsps=%d, h=%.2f, BT=%.2f, span=%d\n', ...
    Rs, Fs, nsps, h, BT, span);
fprintf('Tone spacing=%.0f Hz, Tone LPF fc=%.0f Hz\n', tone_spacing, Fc_tone);
fprintf('Delays: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n', ...
    delay_gauss, delay_ch, delay_tone, total_delay);
fprintf('N_pre=%d, N_post=%d, Nsym_valid=%d, First sample=%d, Last sample=%d\n\n', ...
    N_pre, N_post, Nsym, sample_idx(1), sample_idx(end));

%% ========================================================================
% 2. Gray encoding/decoding mapping
%% ========================================================================
gray_enc = [0; 1; 3; 2];   % 00->0, 01->1, 11->3, 10->2
gry2nat = zeros(4,1);
for i = 0:3
    g = gray_enc(i+1);
    if g <= 3 && g >= 0 && g == round(g)
        gry2nat(g+1) = i;
    end
end
for i = 0:3
    if gry2nat(gray_enc(i+1)+1) ~= i
        error('Gray 映射不一致！');
    end
end

% Frequency indices
freq_no = [-3; -1; 1; 3];
tone_freq = freq_no * h * Rs / 2;

%% ========================================================================
% 3. Helper function: generate GFSK signal
%% ========================================================================
    function s = generate_gfsk(sym_seq)
        Nsym_in = length(sym_seq);
        Ns_in = Nsym_in * nsps;
        
        sym_gray_in = gray_enc(sym_seq + 1);
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
    end

%% ========================================================================
% 4. Helper function: Tone-Mixer soft-decision metric acquisition
%% ========================================================================
    function branch_metric = detect_tonemixer(r)
        Ns_r = length(r);
        t = (0:Ns_r-1)' / Fs;
        branch_metric = zeros(M, Nsym);
        for m = 1:M
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            y_lpf = filter(tone_coeffs, 1, y_mix);
            y_sample = y_lpf(sample_idx);
            branch_metric(m, :) = abs(y_sample).';
        end
    end

%% ========================================================================
% 5. Helper function: 4-state Viterbi soft-decision sequence detection
%% ========================================================================
    function det_gray = viterbi_decode(obs_matrix)
        % obs_matrix: M x Nsym, Per symbolEachBranch'sMatchDegree（Larger means moreMatch）
        [M_v, T] = size(obs_matrix);
        
        % Path metric (pm) WithTracebackMatrix (back)
        pm = zeros(M_v, T);
        back = zeros(M_v, T);
        
        % Initialize：t=1
        pm(:, 1) = obs_matrix(:, 1);
        
        % Forward recursion：t=2..T
        % State = WhenPrevious symbol value；AllAll transitions allowed（4GFSKUnconstrained）
        for t = 2:T
            for s = 1:M_v
                % FromAllPredecessorState transitionToWhenPrevious state s cumulative metric
                cand = pm(:, t-1) + obs_matrix(s, t);
                [pm(s, t), back(s, t)] = max(cand);
            end
            % NormalizationPrevent numerical overflow（Does not affectOptimal path）
            pm(:, t) = pm(:, t) - max(pm(:, t));
        end
        
        % Traceback
        det_gray = zeros(T, 1);
        [~, det_gray(end)] = max(pm(:, end));
        for t = T-1:-1:1
            det_gray(t) = back(det_gray(t+1), t+1);
        end
        det_gray = det_gray - 1;  % ConvertAs 0-based
    end

%% ========================================================================
% 6. Main simulation: Eb/N0 sweep (hard decision vs Viterbi)
%% ========================================================================
EbN0_lin = 10.^(EbN0_dB/10);
BER_hard = zeros(size(EbN0_dB));
BER_vit  = zeros(size(EbN0_dB));
SER_hard = zeros(size(EbN0_dB));
SER_vit  = zeros(size(EbN0_dB));

tic;
for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    N0 = 8 / ebno;           % Eb = 8 (nsps=16, k=2)
    noise_var = N0;          % Complex noise total variance
    
    bit_err_hard = 0; bit_err_vit = 0;
    sym_err_hard = 0; sym_err_vit = 0;
    
    for sim = 1:Nsim
        rng(idx*100 + sim);
        
        % GenerateRandomSymbol
        sym_tx = randi([0, M-1], Nsym_total, 1);
        s = generate_gfsk(sym_tx);
        
        % Power check
        sig_pow = mean(abs(s).^2);
        if abs(sig_pow - 1) > 0.1
            warning('Signal power deviation: %.3f', sig_pow);
        end
        
        % AWGN + channel filter
        noise = sqrt(noise_var/2) * (randn(Ns_total, 1) + 1j*randn(Ns_total, 1));
        r_ch = filter(ch_filter, s + noise);
        
        % AcquireSoft decisionBranch metric
        branch_metric = detect_tonemixer(r_ch);
        
        % ---- Hard decision：Per-symbol maximum magnitude ----
        [~, det_gray_hard] = max(branch_metric, [], 1);
        det_gray_hard = det_gray_hard(:) - 1;
        det_sym_hard = gry2nat(det_gray_hard + 1);
        
        % ---- Viterbi sequence detection ----
        det_gray_vit = viterbi_decode(branch_metric);
        det_sym_vit = gry2nat(det_gray_vit + 1);
        
        % ValidSymbolComparison
        sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
        
        sym_err_hard = sym_err_hard + sum(det_sym_hard ~= sym_tx_valid);
        sym_err_vit  = sym_err_vit  + sum(det_sym_vit  ~= sym_tx_valid);
        
        % BitError
        for i = 1:Nsym
            bit_err_hard = bit_err_hard + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_hard(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_hard(i), 1));
            bit_err_vit = bit_err_vit + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_vit(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_vit(i), 1));
        end
    end
    
    BER_hard(idx) = bit_err_hard / (Nsym * k * Nsim);
    BER_vit(idx)  = bit_err_vit  / (Nsym * k * Nsim);
    SER_hard(idx) = sym_err_hard / (Nsym * Nsim);
    SER_vit(idx)  = sym_err_vit  / (Nsym * Nsim);
    
    fprintf('Eb/N0=%5.1f dB | Hard BER=%.4e | Vit BER=%.4e | Gain=%.2f dB\n', ...
        EbN0_dB(idx), BER_hard(idx), BER_vit(idx), ...
        10*log10(BER_hard(idx)/BER_vit(idx)));
end
fprintf('Total simulation time: %.2f s\n', toc);

%% ========================================================================
% 7. Theoretical BER reference (M-ary orthogonal FSK coherent detection)
%% ========================================================================
BER_theory = zeros(size(EbN0_dB));
for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    sqrt_term = sqrt(2 * ebno * log2(M));
    y_min = -5*sqrt_term - 10;
    y_max = 10;
    if y_min < -50, y_min = -50; end
    Ny = 2000;
    y = linspace(y_min, y_max, Ny);
    phi_y = (1/sqrt(2*pi)) * exp(-y.^2/2);
    Q_term = qfunc(y + sqrt_term);
    integrand = phi_y .* (1 - Q_term).^(M-1);
    P_s = 1 - trapz(y, integrand);
    BER_theory(idx) = P_s / log2(M);
end

BER_bound = (M-1)/log2(M) * qfunc(sqrt(EbN0_lin * log2(M)));

%% ========================================================================
% 8. Visualization
%% ========================================================================

% Figure 1: BER/SER comparison
figure('Name', 'BER Comparison: Hard vs Viterbi', 'Position', [100 100 900 500]);
semilogy(EbN0_dB, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard Decision');
hold on;
semilogy(EbN0_dB, BER_vit, 'r^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Viterbi Sequence');
semilogy(EbN0_dB, BER_theory, 'g-', 'LineWidth', 1, 'DisplayName', 'Theory (orthog. MFSK)');
semilogy(EbN0_dB, BER_bound, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Union Bound');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
legend('Location', 'southwest');
title(sprintf('4-ary GFSK: Hard Decision vs Viterbi (h=%.1f, BT=%.1f, nsps=%d)', h, BT, nsps));
axis([min(EbN0_dB) max(EbN0_dB) 1e-4 1]);

% Figure 2: Viterbi gain (dB)
figure('Name', 'Viterbi Coding Gain', 'Position', [150 150 600 400]);
gain_dB = 10*log10(BER_hard ./ BER_vit);
gain_dB(BER_vit == 0) = 0;  % Avoid division by zero
plot(EbN0_dB, gain_dB, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Viterbi Gain (dB) = 10*log_{10}(BER_{hard}/BER_{vit})');
title('Viterbi Sequence Detection Gain over Hard Decision');

% Figure 3: Single waveform illustration (noiseless, last EbN0 point)
figure('Name', 'Signal & Viterbi Detection', 'Position', [200 200 1000 700]);

rng(123);
sym_demo = randi([0, M-1], Nsym_total, 1);
s_demo = generate_gfsk(sym_demo);
r_demo = filter(ch_filter, s_demo);
branch_demo = detect_tonemixer(r_demo);

% Hard decision
det_gray_hard_demo = zeros(Nsym, 1);
for t = 1:Nsym
    [~, det_gray_hard_demo(t)] = max(branch_demo(:, t));
end
det_gray_hard_demo = det_gray_hard_demo - 1;
det_sym_hard_demo = gry2nat(det_gray_hard_demo + 1);

% Viterbi
det_gray_vit_demo = viterbi_decode(branch_demo);
det_sym_vit_demo = gry2nat(det_gray_vit_demo + 1);

sym_demo_valid = sym_demo(N_pre+1 : N_pre+Nsym);

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

% 3.2 4-branch Metric
subplot(3,2,2);
plot(t_sym*1e3, branch_demo.', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Branch Metric');
legend('Branch 0', 'Branch 1', 'Branch 2', 'Branch 3', 'Location', 'best');
title('Soft Branch Metrics (Tone-Mixer Output)');
grid on; xlim([0 5]);

% 3.3 Hard decision vs transmitted
subplot(3,2,3);
stem(1:50, sym_demo_valid(1:50), 'b', 'LineWidth', 1.5, 'DisplayName', 'Tx');
hold on;
stem(1:50, det_sym_hard_demo(1:50), 'r--', 'LineWidth', 1, 'DisplayName', 'Hard');
xlabel('Symbol Index'); ylabel('Symbol Value');
title('Hard Decision vs Transmitted (First 50)');
legend('Location', 'best'); grid on;

% 3.4 Viterbi vs Transmit
subplot(3,2,4);
stem(1:50, sym_demo_valid(1:50), 'b', 'LineWidth', 1.5, 'DisplayName', 'Tx');
hold on;
stem(1:50, det_sym_vit_demo(1:50), 'g--', 'LineWidth', 1, 'DisplayName', 'Viterbi');
xlabel('Symbol Index'); ylabel('Symbol Value');
title('Viterbi vs Transmitted (First 50)');
legend('Location', 'best'); grid on;

% 3.5 Hard decisionError locations
subplot(3,2,5);
err_hard = (det_sym_hard_demo ~= sym_demo_valid);
err_vit  = (det_sym_vit_demo  ~= sym_demo_valid);
plot(1:100, err_hard(1:100), 'ro', 'MarkerSize', 4, 'DisplayName', 'Hard');
hold on;
plot(1:100, err_vit(1:100)*1.2, 'g^', 'MarkerSize', 4, 'DisplayName', 'Viterbi');
xlabel('Symbol Index'); ylabel('Error Flag');
title('Error Locations: Hard vs Viterbi (First 100)');
legend('Location', 'best'); grid on; ylim([0 1.5]);

% 3.6 Constellation
subplot(3,2,6);
plot(real(r_demo), imag(r_demo), 'b.', 'MarkerSize', 3);
hold on;
plot(real(r_demo(sample_idx)), imag(r_demo(sample_idx)), 'ro', 'MarkerSize', 6);
xlabel('In-Phase'); ylabel('Quadrature');
title('Constellation at Sampling Points');
axis equal; grid on;

%% ========================================================================
% 9. Results summary
%% ========================================================================
fprintf('\n========== RESULT SUMMARY ==========\n');
fprintf('Eb/N0(dB) | Hard BER  |  Vit BER  | Vit Gain(dB) | Theory BER\n');
fprintf('----------|-----------|-----------|--------------|------------\n');
for idx = 1:length(EbN0_dB)
    gain = 10*log10(BER_hard(idx)/BER_vit(idx));
    fprintf('%7.1f   | %.4e | %.4e |    %6.2f    |  %.4e\n', ...
        EbN0_dB(idx), BER_hard(idx), BER_vit(idx), gain, BER_theory(idx));
end
fprintf('\nSimulation complete.\n');

end
