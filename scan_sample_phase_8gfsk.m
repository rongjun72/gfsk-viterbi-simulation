function scan_sample_phase_8gfsk()
% scan_sample_phase_8gfsk.m
% Scan 8-ary GFSK sampling phase offset, find optimal sampling point
%
% Scan range: delta = -3, -2, -1, 0, +1, +2, +3 samples
% For each delta:
%   1. Recompute ISI reference templates (at offset sampling points)
%   2. Hard-decision BER
%   3. Viterbi BER
%   4. Transmitted branch discrimination ratio (Tx / max-adjacent)

%% ========================================================================
% 0. Parameter configuration (fully consistent with gfsk_8ary_viterbi_isi)
%% ========================================================================
Rs      = 1e3;
Fs      = 16e3;
nsps    = Fs/Rs;
M       = 8;
k       = log2(M);
h       = 1.0;
BT      = 0.5;
span    = 4;
Nsym    = 10000;  % Long enough to statistics BER

% Filter design
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

Fp = 4.5e3;   Fs_stop = 5.5e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;
delay_ch = grpdelay(ch_filter.Coefficients,1,1)+0;

tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(24, Fc_tone/(Fs/2), 'low', chebwin(25, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = N_pre;
Nsym_total = Nsym + N_pre + N_post;
Ns_total   = Nsym_total * nsps;

% Gray encoding/decoding
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7, gry2nat(gray_enc(i+1)+1) = i; end
freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== 8-GFSK Sample Phase Scan ===\n');
fprintf('Parameters: Rs=%d, Fs=%d, nsps=%d, h=%.2f, BT=%.2f, span=%d\n', Rs, Fs, nsps, h, BT, span);
fprintf('Total delay = %d samples = %.2f symbols\n', total_delay, total_delay/nsps);
fprintf('Base sample: (k-1)*%d + %d + %d\n', nsps, nsps/2, total_delay);
fprintf('Scanning delta = -3 ~ +3 samples\n\n');

%% ========================================================================
% 1. Helper functions
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

    function bm = measure_tonemixer(r, idx)
        Ns_r = length(r);
        t = (0:Ns_r-1)' / Fs;
        bm = zeros(M, length(idx));
        for m = 1:M
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            y_lpf = filter(tone_coeffs, 1, y_mix);
            bm(m, :) = abs(y_lpf(idx)).';
        end
    end

    function det_gray = viterbi_decode(obs_matrix, ref_m)
        [~, T] = size(obs_matrix);
        N_s = M;
        pm = zeros(N_s, T);
        back = zeros(N_s, T);
        
        % t=1: preamble prev=0
        for curr_g = 0:M-1
            s_idx = curr_g + 1;
            ref = squeeze(ref_m(0+1, curr_g+1, :));
            obs = obs_matrix(:, 1);
            n_obs = norm(obs); n_ref = norm(ref);
            if n_obs > 1e-6 && n_ref > 1e-6
                obs_n = obs / n_obs; ref_n = ref / n_ref;
                branch = obs_n' * ref_n;
            else
                branch = 0;
            end
            pm(s_idx, 1) = branch;
        end
        
        % t=2:T
        for t = 2:T
            for curr_g = 0:M-1
                s_prime = curr_g + 1;
                best_val = -inf; best_prev = 1;
                for prev_g = 0:M-1
                    s = prev_g + 1;
                    ref = squeeze(ref_m(prev_g+1, curr_g+1, :));
                    obs = obs_matrix(:, t);
                    n_obs = norm(obs); n_ref = norm(ref);
                    if n_obs > 1e-6 && n_ref > 1e-6
                        obs_n = obs / n_obs; ref_n = ref / n_ref;
                        branch = obs_n' * ref_n;
                    else
                        branch = 0;
                    end
                    val = pm(s, t-1) + branch;
                    if val > best_val
                        best_val = val; best_prev = s;
                    end
                end
                pm(s_prime, t) = best_val;
                back(s_prime, t) = best_prev;
            end
            pm(:, t) = pm(:, t) - max(pm(:, t));
        end
        
        det_gray = zeros(T, 1);
        [~, s_final] = max(pm(:, end));
        det_gray(end) = s_final - 1;
        for t = T-1:-1:1
            s_prev = back(s_final, t+1);
            det_gray(t) = s_prev - 1;
            s_final = s_prev;
        end
    end

%% ========================================================================
% 2. Generate fixed noiseless test signal
%% ========================================================================
rng(42);
sym_test = [zeros(N_pre, 1); randi([0, M-1], Nsym, 1); zeros(N_post, 1)];
s_test = generate_gfsk(sym_test);
r_test = filter(ch_coeffs, 1, s_test);

% Verify power
sig_pow = mean(abs(s_test).^2);
fprintf('Signal power = %.4f\n\n', sig_pow);

%% ========================================================================
% 3. Scan sampling phase offset
%% ========================================================================
deltas = -3:3;
N_delta = length(deltas);

BER_hard   = zeros(1, N_delta);
BER_vit    = zeros(1, N_delta);
disc_ratio = zeros(1, N_delta);  % Transmitted branch / second largest branch

% Pre-calculate template (for each delta, but use inner loop)
N_guard = 12;

for di = 1:N_delta
    delta = deltas(di);
    fprintf('--- delta = %+d ---\n', delta);
    
    % Sampling point（After offset）
    sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay + delta;
    
    % MeasureBranch metric
    bm = measure_tonemixer(r_test, sample_idx);
    
    % Hard decision
    [~, det_gray_hard] = max(bm, [], 1);
    det_gray_hard = det_gray_hard(:) - 1;
    det_sym_hard = gry2nat(det_gray_hard + 1);
    
    % Calculate discrimination ratio: transmitted branch / maximum interference branch
    tx_symbols = sym_test(N_pre+1 : N_pre+Nsym);
    tx_gray = gray_enc(tx_symbols + 1);
    ratios = zeros(Nsym, 1);
    for t = 1:Nsym
        tx_val = bm(tx_gray(t)+1, t);
        others = bm(:, t);
        others(tx_gray(t)+1) = -inf;  % Exclude transmitted branch
        max_other = max(others);
        if max_other > 1e-6
            ratios(t) = tx_val / max_other;
        else
            ratios(t) = inf;
        end
    end
    disc_ratio(di) = mean(ratios);
    
    % Re-precomputeCalculateISIReference template（At offset sampling point）
    ref_metric = zeros(M, M, M);
    for prev_g = 0:7
        for curr_g = 0:7
            prev_nat = gry2nat(prev_g + 1);
            curr_nat = gry2nat(curr_g + 1);
            sym_seq = [zeros(N_guard, 1); prev_nat; curr_nat; zeros(N_guard, 1)];
            s = generate_gfsk(sym_seq);
            s_ch = filter(ch_coeffs, 1, s);
            k_curr = N_guard + 2;
            idx_curr = (k_curr - 1) * nsps + nsps/2 + total_delay + delta;
            ref_metric(prev_g+1, curr_g+1, :) = measure_tonemixer(s_ch, idx_curr);
        end
    end
    
    % ViterbiDecode
    det_gray_vit = viterbi_decode(bm, ref_metric);
    det_sym_vit = gry2nat(det_gray_vit + 1);
    
    % BERCalculate
    sym_valid = sym_test(N_pre+1 : N_pre+Nsym);
    BER_hard(di) = sum(det_sym_hard ~= sym_valid) / Nsym;
    BER_vit(di)  = sum(det_sym_vit  ~= sym_valid) / Nsym;
    
    fprintf('  Hard BER = %.4e | Vit BER = %.4e | Disc Ratio = %.3f\n', ...
        BER_hard(di), BER_vit(di), disc_ratio(di));
end

%% ========================================================================
% 4. Results summary and visualization
%% ========================================================================
fprintf('\n========== RESULT SUMMARY ==========\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Delta', 'Hard BER', 'Vit BER', 'Disc Ratio');
fprintf('%s\n', repmat('-', 1, 55));
for di = 1:N_delta
    fprintf('%+8d | %.4e | %.4e | %12.3f\n', ...
        deltas(di), BER_hard(di), BER_vit(di), disc_ratio(di));
end

[min_hard, idx_hard] = min(BER_hard);
[min_vit, idx_vit]   = min(BER_vit);
[max_disc, idx_disc] = max(disc_ratio);

fprintf('\nBest Hard BER:   delta=%+d, BER=%.4e\n', deltas(idx_hard), min_hard);
fprintf('Best Viterbi:    delta=%+d, BER=%.4e\n', deltas(idx_vit), min_vit);
fprintf('Best Disc Ratio: delta=%+d, ratio=%.3f\n', deltas(idx_disc), max_disc);

% Current (delta=0) comparison
idx0 = find(deltas == 0);
fprintf('\nCurrent (delta=0): Hard=%.4e, Vit=%.4e, Disc=%.3f\n', ...
    BER_hard(idx0), BER_vit(idx0), disc_ratio(idx0));

% PlotFigure
figure('Name', 'Sample Phase Scan', 'Position', [100 100 1200 400]);

subplot(1, 3, 1);
plot(deltas, BER_hard*100, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard BER');
hold on;
plot(deltas, BER_vit*100, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Vit BER');
xline(0, 'k--', 'Alpha', 0.5, 'HandleVisibility', 'off');
xlabel('Phase Offset (samples)');
ylabel('BER (%)');
title('BER vs Sampling Phase');
legend('Location', 'best');
grid on;

subplot(1, 3, 2);
semilogy(deltas, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
semilogy(deltas, BER_vit, 'r^-', 'LineWidth', 1.5, 'MarkerSize', 8);
xline(0, 'k--', 'Alpha', 0.5, 'HandleVisibility', 'off');
xlabel('Phase Offset (samples)');
ylabel('BER (log)');
title('BER (log scale)');
grid on;

subplot(1, 3, 3);
plot(deltas, disc_ratio, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
xline(0, 'k--', 'Alpha', 0.5, 'HandleVisibility', 'off');
[max_disc_line, max_idx] = max(disc_ratio);
plot(deltas(max_idx), max_disc_line, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Phase Offset (samples)');
ylabel('Tx / Max-Adjacent Ratio');
title('Discrimination Ratio');
grid on;

fprintf('\nScan complete.\n');
end
