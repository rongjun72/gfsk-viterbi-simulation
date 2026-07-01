function gfsk_4ary_viterbi_64state()
% gfsk_4ary_viterbi_64state.m
% 4-ary GFSK coherent demodulation + 64-state ISI-aware Viterbi sequence detection
% State = previous 3 symbols (prev3, prev2, prev1) Gray encoding, 4^3 = 64 states total
% Transition: (prev3, prev2, prev1) -> (prev2, prev1, curr)
% Branch metric based on 4-symbol (prev3, prev2, prev1, curr) ISI reference template
%
% Architecture:
%   1. Tx: continuous-phase GFSK, gaussdesign(BT,span,sps) Gaussian pulse shaping
%   2. Channel: AWGN + 80dB out-of-band rejection channel filter
%   3. Receiver frontend：4-branch tone-mixer coherent detection + Chebyshev window LPF
%   4. Receiver backend：64-State Viterbi Soft decisionSequence detection（3-symbol memory）
%   5. Comparison: per-symbol hard decision vs 4-state Viterbi vs 64-state Viterbi

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
ch_coeffs = ch_filter.Coefficients;

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

fprintf('=== 4-ary GFSK + 64-State ISI Viterbi ===\n');
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
% 5. Precompute: 64-state ISI reference templates
% (prev3_gray, prev2_gray, prev1_gray, curr_gray, branch_idx)
% curr symbol is at position guard+4 in sequence
%% ========================================================================
fprintf('--- Pre-computing 64-state ISI reference templates ---\n');
N_guard = 12;
ref_metric = zeros(4, 4, 4, 4, 4);  % 256 templates, 4-branch each

for prev3_g = 0:3
    for prev2_g = 0:3
        for prev1_g = 0:3
            for curr_g = 0:3
                prev3_n = gry2nat(prev3_g + 1);
                prev2_n = gry2nat(prev2_g + 1);
                prev1_n = gry2nat(prev1_g + 1);
                curr_n = gry2nat(curr_g + 1);
                
                sym_seq = [zeros(N_guard, 1); prev3_n; prev2_n; prev1_n; curr_n; zeros(N_guard, 1)];
                s = generate_gfsk(sym_seq);
                s_ch = filter(ch_coeffs, 1, s);
                
                k_curr = N_guard + 4;  % curr is the (guard+4)-th symbol
                idx_curr = (k_curr - 1) * nsps + nsps/2 + total_delay;
                
                bm = measure_tonemixer(s_ch, idx_curr);
                ref_metric(prev3_g+1, prev2_g+1, prev1_g+1, curr_g+1, :) = bm;
            end
        end
    end
end

fprintf('Reference templates computed: 256 x 4-branch.\n\n');

%% ========================================================================
% 6. Helper function: 64-state Viterbi soft-decision sequence detection
% State = (prev3, prev2, prev1) in Gray, 4^3 = 64 states
%% ========================================================================
    function det_gray = viterbi_64state(obs_matrix)
        [M_v, T] = size(obs_matrix);
        N_s = 64;  % 4^3 states
        
        pm = -inf(N_s, T);
        back = zeros(N_s, T);
        
        % t=1: preamble = (0,0,0)
        for curr = 0:3
            s_prime = 0*16 + 0*4 + curr + 1;  % (0,0,curr), 1-based
            ref = squeeze(ref_metric(0+1, 0+1, 0+1, curr+1, :));
            obs = obs_matrix(:, 1);
            n_obs = norm(obs); n_ref = norm(ref);
            if n_obs > 1e-6 && n_ref > 1e-6
                pm(s_prime, 1) = (obs/n_obs)' * (ref/n_ref);
            else
                pm(s_prime, 1) = 0;
            end
        end
        
        % t=2:T: forward recursion
        for t = 2:T
            for curr = 0:3
                for prev2 = 0:3
                    for prev1 = 0:3
                        s_prime = prev2*16 + prev1*4 + curr + 1;  % 1-based
                        best_val = -inf;
                        best_prev_state = 1;
                        for prev3 = 0:3
                            s = prev3*16 + prev2*4 + prev1 + 1;  % 1-based
                            if pm(s, t-1) == -inf
                                continue;
                            end
                            ref = squeeze(ref_metric(prev3+1, prev2+1, prev1+1, curr+1, :));
                            obs = obs_matrix(:, t);
                            n_obs = norm(obs); n_ref = norm(ref);
                            if n_obs > 1e-6 && n_ref > 1e-6
                                branch = (obs/n_obs)' * (ref/n_ref);
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
                end
            end
            % Normalize to prevent overflow
            pm(:, t) = pm(:, t) - max(pm(:, t));
        end
        
        % Traceback
        det_gray = zeros(T, 1);
        [~, s_final] = max(pm(:, end));
        for t = T:-1:1
            curr = mod(s_final-1, 4);
            det_gray(t) = curr;
            if t > 1
                s_final = back(s_final, t);
            end
        end
    end

%% ========================================================================
% 6.5 Noiseless Viterbi self-check (64-state vs hard decision vs 4-state)
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

% 64-State Viterbi
det_gray_vit64_test = viterbi_64state(bm_test);
det_sym_vit64_test = gry2nat(det_gray_vit64_test + 1);

sym_test_valid = sym_test(N_pre+1 : N_pre+Nsym);
BER_hard_test = sum(det_sym_hard_test ~= sym_test_valid) / Nsym;
BER_vit64_test = sum(det_sym_vit64_test ~= sym_test_valid) / Nsym;

fprintf('  Hard decision: BER = %.4e\n', BER_hard_test);
fprintf('  64-State Viterbi-ISI: BER = %.4e\n', BER_vit64_test);
if BER_vit64_test > 0
    fprintf('  WARNING: 64-state Viterbi has errors even without noise!\n');
else
    fprintf('  64-state Viterbi passes noiseless check.\n');
end
fprintf('\n');

%% ========================================================================
% 7. Main simulation: Eb/N0 sweep (hard decision vs 64-state Viterbi)
%% ========================================================================
EbN0_lin = 10.^(EbN0_dB/10);
BER_hard = zeros(size(EbN0_dB));
BER_vit64  = zeros(size(EbN0_dB));
SER_hard = zeros(size(EbN0_dB));
SER_vit64  = zeros(size(EbN0_dB));

tic;
for idx = 1:length(EbN0_dB)
    ebno = EbN0_lin(idx);
    N0 = 8 / ebno;           % Eb = 8 (nsps=16, k=2)
    noise_var = N0;          % Complex noise total variance
    
    bit_err_hard = 0; bit_err_vit64 = 0;
    sym_err_hard = 0; sym_err_vit64 = 0;
    
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
        
        % ---- 64-State Viterbi sequence detection ----
        det_gray_vit64 = viterbi_64state(branch_metric);
        det_sym_vit64 = gry2nat(det_gray_vit64 + 1);
        
        % ValidSymbolComparison
        sym_tx_valid = sym_tx(N_pre+1 : N_pre+Nsym);
        
        sym_err_hard = sym_err_hard + sum(det_sym_hard ~= sym_tx_valid);
        sym_err_vit64  = sym_err_vit64  + sum(det_sym_vit64  ~= sym_tx_valid);
        
        % BitError
        for i = 1:Nsym
            bit_err_hard = bit_err_hard + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_hard(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_hard(i), 1));
            bit_err_vit64 = bit_err_vit64 + ...
                (bitget(sym_tx_valid(i), 2) ~= bitget(det_sym_vit64(i), 2)) + ...
                (bitget(sym_tx_valid(i), 1) ~= bitget(det_sym_vit64(i), 1));
        end
    end
    
    BER_hard(idx) = bit_err_hard / (Nsym * k * Nsim);
    BER_vit64(idx)  = bit_err_vit64  / (Nsym * k * Nsim);
    SER_hard(idx) = sym_err_hard / (Nsym * Nsim);
    SER_vit64(idx)  = sym_err_vit64  / (Nsym * Nsim);
    
    gain = 10*log10(BER_hard(idx)/BER_vit64(idx));
    fprintf('Eb/N0=%5.1f dB | Hard BER=%.4e | 64-Vit BER=%.4e | Gain=%.2f dB\n', ...
        EbN0_dB(idx), BER_hard(idx), BER_vit64(idx), gain);
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
figure('Name', 'BER: Hard vs 64-State Viterbi', 'Position', [100 100 900 500]);
semilogy(EbN0_dB, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard Decision');
hold on;
semilogy(EbN0_dB, BER_vit64, 'r^-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '64-State Viterbi-ISI');
semilogy(EbN0_dB, BER_theory, 'g-', 'LineWidth', 1, 'DisplayName', 'Theory (orthog. MFSK)');
semilogy(EbN0_dB, BER_bound, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Union Bound');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate');
legend('Location', 'southwest');
title(sprintf('4-ary GFSK: Hard vs 64-State Viterbi-ISI (h=%.1f, BT=%.1f, nsps=%d)', h, BT, nsps));
axis([min(EbN0_dB) max(EbN0_dB) 1e-4 1]);

% Figure 2: 64-State Viterbi gain (dB)
figure('Name', '64-State Viterbi Gain', 'Position', [150 150 600 400]);
gain_dB = 10*log10(BER_hard ./ BER_vit64);
gain_dB(BER_vit64 == 0) = 0;
plot(EbN0_dB, gain_dB, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('64-State Viterbi Gain (dB)');
title('64-State Viterbi-ISI Gain over Hard Decision');

% Figure 3: Single waveform illustration (noiseless, last EbN0 point)
figure('Name', 'Signal & 64-State Viterbi Detection', 'Position', [200 200 1000 700]);

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

% 64-State Viterbi
det_gray_vit64_demo = viterbi_64state(branch_demo);
det_sym_vit64_demo = gry2nat(det_gray_vit64_demo + 1);

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

% 3.2 4-branch soft metrics
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

% 3.4 64-State Viterbi vs transmitted
subplot(3,2,4);
stem(1:50, sym_demo_valid(1:50), 'b', 'LineWidth', 1.5, 'DisplayName', 'Tx');
hold on;
stem(1:50, det_sym_vit64_demo(1:50), 'g--', 'LineWidth', 1, 'DisplayName', '64-Viterbi');
xlabel('Symbol Index'); ylabel('Symbol Value');
title('64-State Viterbi vs Transmitted (First 50)');
legend('Location', 'best'); grid on;

% 3.5 Error location comparison
subplot(3,2,5);
err_hard = (det_sym_hard_demo ~= sym_demo_valid);
err_vit64 = (det_sym_vit64_demo ~= sym_demo_valid);
plot(1:100, err_hard(1:100), 'ro', 'MarkerSize', 4, 'DisplayName', 'Hard');
hold on;
plot(1:100, err_vit64(1:100)*1.2, 'g^', 'MarkerSize', 4, 'DisplayName', '64-Viterbi');
xlabel('Symbol Index'); ylabel('Error Flag');
title('Error Locations: Hard vs 64-State Viterbi (First 100)');
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
fprintf('Eb/N0(dB) | Hard BER  | 64-Vit BER | 64-Vit Gain(dB) | Theory BER\n');
fprintf('----------|-----------|------------|-----------------|------------\n');
for idx = 1:length(EbN0_dB)
    gain = 10*log10(BER_hard(idx)/BER_vit64(idx));
    fprintf('%7.1f   | %.4e |  %.4e  |     %6.2f     |  %.4e\n', ...
        EbN0_dB(idx), BER_hard(idx), BER_vit64(idx), gain, BER_theory(idx));
end
fprintf('\nSimulation complete.\n');

end
