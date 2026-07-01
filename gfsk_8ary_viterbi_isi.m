function gfsk_8ary_viterbi_isi()
% gfsk_8ary_viterbi_isi.m
% 8-ary GFSK coherent demodulation + ISI-aware Viterbi sequence detection
%
% Architecture:
%   1. Tx: continuous-phase 8-GFSK, gaussdesign(BT,span,sps) Gaussian pulse shaping
%   2. Channel: AWGN + 80dB out-of-band rejection channel filter
%   3. Rx front-end: 8-branch tone-mixer coherent detection + Chebyshev window LPF
%   4. Rx back-end: 8-state ISI-aware Viterbi sequence detection (1-symbol memory)
%
% ISI-aware design:
%   - Precompute 64 standard noiseless 8-GFSK waveforms for (prev_gray, curr_gray) combinations
%   - Measure 8-branch tone-mixer output at current symbol midpoint as reference template
%   - Normalize reference templates, Viterbi branch metric = inner product of observation and reference (cosine similarity)
%   - State = curr Gray code; transition prev→curr
%   - Normalization prevents overflow, full-frame traceback

%% ========================================================================
% 0. Configurable parameters
%% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 8;            % 8-ary
k       = log2(M);      % 3 bits/symbol
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
% Optimization: 8-GFSK outermost tone at 3500 Hz, 99% signal power within ±4000 Hz
Fp = 4.5e3;   Fs_stop = 5.5e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
delay_ch = grpdelay(ch_filter.Coefficients,1,1)+0;
ch_coeffs = ch_filter.Coefficients;

% 1.3 Tone mixer lowpass filter: Chebyshev window, 24-tap, fc=0.75*tone spacing
tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(24, Fc_tone/(Fs/2), 'low', chebwin(25, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% Total delay (for main simulation steady-state calculation)
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

fprintf('=== 8-ary GFSK + ISI-Aware Viterbi Sequence Detection ===\n');
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
% 8-ary 3-bit Gray encoding
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];  % 000->0, 001->1, 010->3, 011->2, 100->6, 101->7, 110->5, 111->4
gry2nat = zeros(8,1);
for i = 0:7
    g = gray_enc(i+1);
    if g <= 7 && g >= 0 && g == round(g)
        gry2nat(g+1) = i;
    end
end
for i = 0:7
    if gry2nat(gray_enc(i+1)+1) ~= i
        error('Gray 映射不一致！');
    end
end

% Frequency indices
freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
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
% 4. Helper function: Tone-Mixer soft-decision metric acquisition (any signal)
%% ========================================================================
    function branch_metric = measure_tonemixer(r, idx)
        Ns_r = length(r);
        t = (0:Ns_r-1)' / Fs;
        branch_metric = zeros(M, length(idx));
        for m = 1:M
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            y_lpf = filter(tone_coeffs, 1, y_mix);
            branch_metric(m, :) = abs(y_lpf(idx)).';
        end
    end

%% ========================================================================
% 5. Precompute: ISI reference templates (prev_gray, curr_gray) 8-branch metrics
% curr symbol is at position guard+2 in sequence (prev, curr)
%% ========================================================================
fprintf('--- Pre-computing 64 ISI reference templates (2-symbol, with ch_filter) ---\n');
N_guard = 12;
ref_metric = zeros(M, M, M);  % (prev_gray, curr_gray, branch_idx)

for prev_g = 0:7
    for curr_g = 0:7
        % ref_metric indexed by GRAY code (consistent with Viterbi states)
        % generate_gfsk expects natural binary, so convert Gray->nat
        prev_nat = gry2nat(prev_g + 1);
        curr_nat = gry2nat(curr_g + 1);
        sym_seq = [zeros(N_guard, 1); prev_nat; curr_nat; zeros(N_guard, 1)];
        s = generate_gfsk(sym_seq);
        s_ch = filter(ch_coeffs, 1, s);  % with ch_filter, consistent with main sim
        
        k_curr = N_guard + 2;  % curr is the (guard+2)-th symbol
        % Sampling point andMain simulationFullyConsistent：Include total_delay
        idx_curr = (k_curr - 1) * nsps + nsps/2 + total_delay;
        
        bm = measure_tonemixer(s_ch, idx_curr);
        ref_metric(prev_g+1, curr_g+1, :) = bm;
    end
end

% Verify reference templates
fprintf('Reference templates computed (64 templates).\n');
count = 0;
for prev_g = 0:7
    for curr_g = 0:7
        ref = squeeze(ref_metric(prev_g+1, curr_g+1, :));
        [~, peak] = max(ref);
        if count < 10  % Only print first 10
            fprintf('  (prev=%d, curr=%d): peak_branch=%d, norm=%.3f\n', ...
                prev_g, curr_g, peak-1, norm(ref));
            count = count + 1;
        end
    end
end
fprintf('  ... (%d more templates)\n\n', 64 - count);

%% ========================================================================
% 6. Helper function: ISI-aware 8-state Viterbi soft-decision sequence detection (1-symbol memory)
% State = curr_gray, 8 states total
% Transition: from state prev to state curr
% Branch metric based on 2-symbol reference template (prev, curr)
%% ========================================================================
    function det_gray = viterbi_decode_isi(obs_matrix)
        [~, T] = size(obs_matrix);  % M=8 branches, T=Nsym symbols
        N_s = M;  % 8 states: curr indexed as curr + 1
        
        pm = zeros(N_s, T);
        back = zeros(N_s, T);
        
        % t=1: PreambleSymbol is 0
        % State s = curr_1，BranchBased on (prev=0, curr_1)
        for curr_g = 0:M-1
            s_idx = curr_g + 1;  % state = curr_g
            ref = squeeze(ref_metric(0+1, curr_g+1, :));
            obs = obs_matrix(:, 1);
            
            n_obs = norm(obs);
            n_ref = norm(ref);
            if n_obs > 1e-6 && n_ref > 1e-6
                obs_n = obs / n_obs;
                ref_n = ref / n_ref;
                branch = obs_n' * ref_n;  % cosine similarity in [0, 1]
            else
                branch = 0;
            end
            
            pm(s_idx, 1) = branch;
        end
        
        % t=2:T: forward recursion
        % New state s' = curr
        % Predecessor state s = prev for all prev
        for t = 2:T
            for curr_g = 0:M-1
                s_prime = curr_g + 1;  % new state = curr
                best_val = -inf;
                best_prev_state = 1;
                
                for prev_g = 0:M-1
                    s = prev_g + 1;  % prev state = prev
                    
                    ref = squeeze(ref_metric(prev_g+1, curr_g+1, :));
                    obs = obs_matrix(:, t);
                    
                    n_obs = norm(obs);
                    n_ref = norm(ref);
                    if n_obs > 1e-6 && n_ref > 1e-6
                        obs_n = obs / n_obs;
                        ref_n = ref / n_ref;
                        branch = obs_n' * ref_n;  % cosine similarity [0,1]
                    else
                        branch = 0;
                    end
                    
                    val = pm(s, t-1) + branch;
                    if val > best_val
                        best_val = val;
                        best_prev_state = s;
                    end
                end
                
                pm(s_prime, t) = best_val;
                back(s_prime, t) = best_prev_state;
            end
            % Normalize to prevent numerical overflow
            pm(:, t) = pm(:, t) - max(pm(:, t));
        end
        
        % Traceback: state = curr directly
        det_gray = zeros(T, 1);
        [~, s_final] = max(pm(:, end));
        det_gray(end) = s_final - 1;  % curr component of last state
        
        for t = T-1:-1:1
            s_prev = back(s_final, t+1);
            det_gray(t) = s_prev - 1;  % curr component of prev state
            s_final = s_prev;
        end
    end

%% ========================================================================
% 6.5 Noiseless Viterbi self-check
%% ========================================================================
fprintf('--- Noiseless Viterbi self-check ---\n');
rng(42);
sym_test = [zeros(N_pre, 1); randi([0, M-1], Nsym, 1); zeros(N_post, 1)];
s_test = generate_gfsk(sym_test);
r_test = filter(ch_coeffs, 1, s_test);
bm_test = measure_tonemixer(r_test, sample_idx);

% Hard decision
det_gray_hard_test = zeros(Nsym, 1);
for t = 1:Nsym
    [~, det_gray_hard_test(t)] = max(bm_test(:, t));
end
det_gray_hard_test = det_gray_hard_test(:) - 1;
det_sym_hard_test = gry2nat(det_gray_hard_test + 1);

% Viterbi
det_gray_vit_test = viterbi_decode_isi(bm_test);
det_sym_vit_test = gry2nat(det_gray_vit_test + 1);

sym_test_valid = sym_test(N_pre+1 : N_pre+Nsym);
BER_hard_test = sum(det_sym_hard_test ~= sym_test_valid) / Nsym;
BER_vit_test = sum(det_sym_vit_test ~= sym_test_valid) / Nsym;

fprintf('  Hard decision: BER = %.4e\n', BER_hard_test);
fprintf('  Viterbi-ISI:   BER = %.4e\n', BER_vit_test);
if BER_vit_test > 0
    fprintf('  WARNING: Viterbi has errors even without noise!\n');
    fprintf('  First 10 errors: Tx=%s, Viterbi=%s\n', ...
        mat2str(sym_test_valid(1:10)'), mat2str(det_sym_vit_test(1:10)'));
else
    fprintf('  Viterbi passes noiseless check.\n');
end
fprintf('\n');

%% ========================================================================
% 7. Main simulation: Eb/N0 sweep (hard decision vs ISI-aware Viterbi)
%% ========================================================================
EbN0_lin = 10.^(EbN0_dB/10);
BER_hard = zeros(size(EbN0_dB));
BER_vit  = zeros(size(EbN0_dB));
SER_hard = zeros(size(EbN0_dB));
SER_vit  = zeros(size(EbN0_dB));

tic;
for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    N0 = (nsps / k) / ebno;  % Eb = nsps/k = 16/3
    noise_var = N0;          % Complex noise total variance
    
    bit_err_hard = 0; bit_err_vit = 0;
    sym_err_hard = 0; sym_err_vit = 0;
    
    for sim = 1:Nsim
        rng(idx*100 + sim);
        
        % GenerateRandomSymbol（Preamble/PostambleFixedAs0）
        sym_tx = [zeros(N_pre, 1); randi([0, M-1], Nsym, 1); zeros(N_post, 1)];
        s = generate_gfsk(sym_tx);
        
        % Power check
        sig_pow = mean(abs(s).^2);
        if abs(sig_pow - 1) > 0.1
            warning('Signal power deviation: %.3f', sig_pow);
        end
        
        % AWGN + channel filter
        noise = sqrt(noise_var/2) * (randn(Ns_total, 1) + 1j*randn(Ns_total, 1));
        r_ch = filter(ch_coeffs, 1, s + noise);
        
        % AcquireSoft decisionBranch metric
        branch_metric = measure_tonemixer(r_ch, sample_idx);
        
        % ---- Hard decision：Per-symbol maximum magnitude ----
        [~, det_gray_hard] = max(branch_metric, [], 1);
        det_gray_hard = det_gray_hard(:) - 1;
        det_sym_hard = gry2nat(det_gray_hard + 1);
        
        % ---- ISI-aware Viterbi sequence detection ----
        det_gray_vit = viterbi_decode_isi(branch_metric);
        det_sym_vit = gry2nat(det_gray_vit + 1);
        
        % ValidSymbolComparison
        sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
        
        sym_err_hard = sym_err_hard + sum(det_sym_hard ~= sym_tx_valid);
        sym_err_vit  = sym_err_vit  + sum(det_sym_vit  ~= sym_tx_valid);
        
        % BitError
        for i = 1:Nsym
            bit_err_hard = bit_err_hard + ...
                (bitget(sym_tx_valid(i), 3) ~= bitget(det_sym_hard(i), 3)) + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_hard(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_hard(i), 1));
            bit_err_vit = bit_err_vit + ...
                (bitget(sym_tx_valid(i), 3) ~= bitget(det_sym_vit(i), 3)) + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_vit(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_vit(i), 1));
        end
    end
    
    BER_hard(idx) = bit_err_hard / (Nsym * k * Nsim);
    BER_vit(idx)  = bit_err_vit  / (Nsym * k * Nsim);
    SER_hard(idx) = sym_err_hard / (Nsym * Nsim);
    SER_vit(idx)  = sym_err_vit  / (Nsym * Nsim);
    
    gain = 10*log10(BER_hard(idx)/BER_vit(idx));
    fprintf('Eb/N0=%5.1f dB | Hard BER=%.4e | Vit BER=%.4e | Gain=%.2f dB\n', ...
        EbN0_dB(idx), BER_hard(idx), BER_vit(idx), gain);
end
fprintf('Total simulation time: %.2f s\n', toc);

%% ========================================================================
% 8. Theoretical BER reference (M-ary orthogonal FSK coherent detection)
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
% 9. Visualization
%% ========================================================================

% Figure 1: BER/SER comparison
figure('Name', 'BER: Hard vs ISI-Aware Viterbi', 'Position', [100 100 900 500]);
semilogy(EbN0_dB, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard Decision');
hold on;
semilogy(EbN0_dB, BER_vit, 'r^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'ISI-Aware Viterbi');
semilogy(EbN0_dB, BER_theory, 'g-', 'LineWidth', 1, 'DisplayName', 'Theory (orthog. MFSK)');
semilogy(EbN0_dB, BER_bound, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Union Bound');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
legend('Location', 'southwest');
title(sprintf('8-ary GFSK: Hard vs ISI-Aware Viterbi (h=%.1f, BT=%.1f, nsps=%d)', h, BT, nsps));
axis([min(EbN0_dB) max(EbN0_dB) 1e-4 1]);

% Figure 2: Viterbi gain (dB)
figure('Name', 'ISI-Aware Viterbi Gain', 'Position', [150 150 600 400]);
gain_dB = 10*log10(BER_hard ./ BER_vit);
gain_dB(BER_vit == 0) = 0;
plot(EbN0_dB, gain_dB, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Viterbi Gain (dB)');
title('ISI-Aware Viterbi Gain over Hard Decision');

% Figure 3: Single waveform illustration (noiseless, last EbN0 point)
figure('Name', 'Signal & ISI-Aware Viterbi Detection', 'Position', [200 200 1000 700]);

rng(123);
sym_demo = [zeros(N_pre, 1); randi([0, M-1], Nsym, 1); zeros(N_post, 1)];
s_demo = generate_gfsk(sym_demo);
r_demo = filter(ch_coeffs, 1, s_demo);
branch_demo = measure_tonemixer(r_demo, sample_idx);

% Hard decision
det_gray_hard_demo = zeros(Nsym, 1);
for t = 1:Nsym
    [~, det_gray_hard_demo(t)] = max(branch_demo(:, t));
end
det_gray_hard_demo = det_gray_hard_demo(:) - 1;
det_sym_hard_demo = gry2nat(det_gray_hard_demo + 1);

% Viterbi
det_gray_vit_demo = viterbi_decode_isi(branch_demo);
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
title('Transmitted 8-GFSK Spectrum');
grid on; xlim([-6 6]);

% 3.2 8-branch Soft metric
subplot(3,2,2);
plot(t_sym*1e3, branch_demo.', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('Branch Metric');
legend('B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'Location', 'best');
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

% 3.4 ISI-aware Viterbi vs Transmit
subplot(3,2,4);
stem(1:50, sym_demo_valid(1:50), 'b', 'LineWidth', 1.5, 'DisplayName', 'Tx');
hold on;
stem(1:50, det_sym_vit_demo(1:50), 'g--', 'LineWidth', 1, 'DisplayName', 'Viterbi-ISI');
xlabel('Symbol Index'); ylabel('Symbol Value');
title('ISI-Aware Viterbi vs Transmitted (First 50)');
legend('Location', 'best'); grid on;

% 3.5 Error location comparison
subplot(3,2,5);
err_hard = (det_sym_hard_demo ~= sym_demo_valid);
err_vit  = (det_sym_vit_demo  ~= sym_demo_valid);
plot(1:100, err_hard(1:100), 'ro', 'MarkerSize', 4, 'DisplayName', 'Hard');
hold on;
plot(1:100, err_vit(1:100)*1.2, 'g^', 'MarkerSize', 4, 'DisplayName', 'Viterbi-ISI');
xlabel('Symbol Index'); ylabel('Error Flag');
title('Error Locations: Hard vs ISI-Aware Viterbi (First 100)');
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
% 10. Results summary
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
