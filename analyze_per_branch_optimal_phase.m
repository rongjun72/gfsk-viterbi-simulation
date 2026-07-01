function analyze_per_branch_optimal_phase()
% analyze_per_branch_optimal_phase.m
% Test: per-branch optimal sampling phase vs global uniform sampling phase
%
% Principle: tones at different frequencies experience different group delays after channel filter + LPF,
%       so optimal sampling phases may differ per branch.
% This script:
%   1. For each branch m, scan all 16 sampling phases within symbol window
%   2. Find peak position for each branch (as the branch's "optimal phase")
%   3. Build "per-branch optimal" metric matrix: branch_m sampled at its optimal phase
%   4. Compare: global uniform sampling vs per-branch optimal sampling BER/discrimination

%% ========================================================================
% 0. Parameter configuration (consistent with main simulation)
%% ========================================================================
Rs      = 1e3;          Fs      = 16e3;         nsps    = Fs/Rs;
M       = 8;            k       = log2(M);      h       = 1.0;
BT      = 0.5;          span    = 4;

% Filters
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

Fp = 4.5e3;   Fs_stop = 5.5e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;   delay_ch = grpdelay(ch_coeffs,1,1)+0;

tone_spacing = h * Rs;   Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(24, Fc_tone/(Fs/2), 'low', chebwin(25, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = N_pre;

% Gray encoding/decoding
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7, gry2nat(gray_enc(i+1)+1) = i; end
freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== Per-Branch Optimal Phase Analysis ===\n');
fprintf('Total delay = %d samples = %.2f symbols\n\n', total_delay, total_delay/nsps);

%% ========================================================================
% 1. Helper functions
%% ========================================================================
    function s = generate_gfsk(sym_seq)
        Nsym_in = length(sym_seq);  Ns_in = Nsym_in * nsps;
        sym_gray_in = gray_enc(sym_seq + 1);
        f_seq = freq_no(sym_gray_in + 1);
        f_up = repelem(f_seq, nsps);
        f_smooth = filter(gauss_filt, 1, f_up);
        dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
        phase = cumsum(dphi);
        s = exp(1j * phase);
        if length(s) < Ns_in, s = [s; zeros(Ns_in - length(s), 1)];
        else, s = s(1:Ns_in); end
    end

    function bm = measure_tonemixer(r, idx)
        Ns_r = length(r);  t = (0:Ns_r-1)' / Fs;
        bm = zeros(M, length(idx));
        for m = 1:M
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            y_lpf = filter(tone_coeffs, 1, y_mix);
            bm(m, :) = abs(y_lpf(idx)).';
        end
    end

%% ========================================================================
% 2. Generate test signals (long sequence for BER stats + short sequence for visualization)
%% ========================================================================
% 2.1 Long sequence (1000 symbols, BER statistics)
rng(42);
Nsym_stat = 1000;
sym_long = [zeros(N_pre,1); randi([0,M-1],Nsym_stat,1); zeros(N_post,1)];
s_long = generate_gfsk(sym_long);
r_long = filter(ch_coeffs, 1, s_long);

% 2.2 Short sequence (24 symbols, per-symbol visualization)
sym_short = [zeros(N_pre,1); [0;1;2;3;4;5;6;7; 7;0;7;0; 0;7;1;6;2;5; 3;3;3;4;4;4]; zeros(N_post,1)];
Nsym_short = 24;
s_short = generate_gfsk(sym_short);
r_short = filter(ch_coeffs, 1, s_short);

fprintf('Generated signals: long=%d symbols, short=%d symbols\n', Nsym_stat, Nsym_short);

%% ========================================================================
% 3. Full-time branch metrics (all samples, no decimation)
%% ========================================================================
t_all = (0:length(r_long)-1)' / Fs;
branch_metric_full_long = zeros(M, length(r_long));
for m = 1:M
    y_mix = r_long .* exp(-1j * 2*pi * tone_freq(m) * t_all);
    y_lpf = filter(tone_coeffs, 1, y_mix);
    branch_metric_full_long(m, :) = abs(y_lpf).';
end

branch_metric_full_short = zeros(M, length(r_short));
t_all_short = (0:length(r_short)-1)' / Fs;
for m = 1:M
    y_mix = r_short .* exp(-1j * 2*pi * tone_freq(m) * t_all_short);
    y_lpf = filter(tone_coeffs, 1, y_mix);
    branch_metric_full_short(m, :) = abs(y_lpf).';
end

%% ========================================================================
% 4. Find "optimal phase" for each branch within symbol window
%% ========================================================================
% Method: for short sequence, analyze 16 samples per branch per symbol
% Find which relative phase (1~16) has the most peaks for each branch

valid_short = N_pre+1 : N_pre+Nsym_short;
branch_peak_pos = zeros(M, Nsym_short);  % Peak value location of each branch at each symbol

for sym_idx = 1:Nsym_short
    k = valid_short(sym_idx);
    % Delay compensation window (16 samples)
    win_start = (k-1)*nsps + total_delay + 1;
    win_end = (k-1)*nsps + total_delay + nsps;
    win = win_start:win_end;
    metrics = branch_metric_full_short(:, win);  % M x 16
    
    for m = 1:M
        [~, peak_idx] = max(metrics(m, :));
        branch_peak_pos(m, sym_idx) = peak_idx;  % 1~16
    end
end

% Statistics for each branch's optimal phase (mode/mean)
fprintf('\n--- Per-Branch Optimal Phase (relative to delay-compensated window) ---\n');
fprintf('%-8s | %-10s | %-10s | %-10s | %-10s\n', 'Branch', 'Mean Pos', 'Median', 'Mode', 'Std');
fprintf('%s\n', repmat('-', 1, 60));

optimal_phase = zeros(M, 1);
for m = 1:M
    mean_pos = mean(branch_peak_pos(m, :));
    median_pos = median(branch_peak_pos(m, :));
    mode_pos = mode(branch_peak_pos(m, :));
    std_pos = std(branch_peak_pos(m, :));
    optimal_phase(m) = round(mean_pos);  % Using mean value as "optimal phase"
    fprintf('%8d | %10.2f | %10.1f | %10d | %10.2f\n', m-1, mean_pos, median_pos, mode_pos, std_pos);
end

%% ========================================================================
% 5. Comparison of three sampling strategies
%% ========================================================================
% Strategy A: global uniform delta=0 (main simulation default)
% Strategy B: global uniform delta=-1 (previously scanned optimal)
% Strategy C: per-branch uses its own optimal phase (theoretical upper bound)

valid_long = N_pre+1 : N_pre+Nsym_stat;
strategies = {'Global delta=0', 'Global delta=-1', 'Per-Branch Optimal'};

% Global sampling point definitions
sample_idx_0 = (valid_long-1)*nsps + nsps/2 + total_delay;       % delta=0
sample_idx_m1 = (valid_long-1)*nsps + nsps/2 + total_delay - 1;  % delta=-1

% Strategy A: global delta=0
bm_A = branch_metric_full_long(:, sample_idx_0);
% Strategy B: global delta=-1
bm_B = branch_metric_full_long(:, sample_idx_m1);

% Strategy C: per-branch optimal phase
bm_C = zeros(M, Nsym_stat);
for m = 1:M
    % Branchmoptimal phase = optimal_phase(m)（At16Position within point window）
    % ConvertAs globalSampling index: Window start + optimal_phase - 1
    sample_idx_per_branch = (valid_long-1)*nsps + total_delay + optimal_phase(m);
    bm_C(m, :) = branch_metric_full_long(m, sample_idx_per_branch);
end

% Hard-decision BER
sym_tx = sym_long(valid_long);
BER_hard_A = 0; BER_hard_B = 0; BER_hard_C = 0;

for t = 1:Nsym_stat
    [~, det_A] = max(bm_A(:, t));   det_A = det_A - 1;   sym_A = gry2nat(det_A + 1);
    [~, det_B] = max(bm_B(:, t));   det_B = det_B - 1;   sym_B = gry2nat(det_B + 1);
    [~, det_C] = max(bm_C(:, t));   det_C = det_C - 1;   sym_C = gry2nat(det_C + 1);
    
    if sym_A ~= sym_tx(t), BER_hard_A = BER_hard_A + 1; end
    if sym_B ~= sym_tx(t), BER_hard_B = BER_hard_B + 1; end
    if sym_C ~= sym_tx(t), BER_hard_C = BER_hard_C + 1; end
end

BER_hard_A = BER_hard_A / Nsym_stat;
BER_hard_B = BER_hard_B / Nsym_stat;
BER_hard_C = BER_hard_C / Nsym_stat;

fprintf('\n--- Hard Decision BER (no noise) ---\n');
SER_vals = [BER_hard_A, BER_hard_B, BER_hard_C];
for i = 1:3
    fprintf('%-20s: SER = %.4e (%.2f%%)\n', strategies{i}, ...
        SER_vals(i), SER_vals(i)*100);
end

%% ========================================================================
% 6. Discrimination ratio comparison
%% ========================================================================
disc_A = compute_disc_ratio(bm_A, sym_tx, gray_enc, gry2nat);
disc_B = compute_disc_ratio(bm_B, sym_tx, gray_enc, gry2nat);
disc_C = compute_disc_ratio(bm_C, sym_tx, gray_enc, gry2nat);

fprintf('\n--- Discrimination Ratio (Tx / Max-Adjacent) ---\n');
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Global delta=0', disc_A.mean, disc_A.std, disc_A.min);
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Global delta=-1', disc_B.mean, disc_B.std, disc_B.min);
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Per-Branch Optimal', disc_C.mean, disc_C.std, disc_C.min);

%% ========================================================================
% 7. Visualization: per-symbol 8-branch time-domain curves + mark different sampling strategies
%% ========================================================================
figure('Name', 'Per-Branch Optimal Phase', 'Position', [50 50 1600 1000]);

n_show = min(8, Nsym_short);  % Show first 8 symbols
for sym_idx = 1:n_show
    k = valid_short(sym_idx);
    win_start = (k-1)*nsps + total_delay + 1;
    win_end = (k-1)*nsps + total_delay + nsps;
    win = win_start:win_end;
    metrics = branch_metric_full_short(:, win);  % 8 x 16
    
    tx_gray = gray_enc(sym_short(k)+1);
    tx_nat = sym_short(k);
    
    subplot(4, 2, sym_idx);
    hold on;
    
    for m = 1:M
        branch_gray = m - 1;
        if branch_gray == tx_gray
            plot(1:nsps, metrics(m,:), 'LineWidth', 2.5, 'Color', 'r', 'DisplayName', sprintf('Tx(B%d)', branch_gray));
        elseif abs(branch_gray - tx_gray) == 1 || abs(branch_gray - tx_gray) == 7
            plot(1:nsps, metrics(m,:), 'LineWidth', 1.2, 'Color', [1 0.6 0], 'DisplayName', sprintf('Adj(B%d)', branch_gray));
        else
            plot(1:nsps, metrics(m,:), 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
    end
    
    % Mark positions of three sampling strategies
    xline(nsps/2, 'b--', 'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    xline(nsps/2 - 1, 'm:', 'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    
    % Mark optimal phase of each branch
    for m = 1:M
        opt_p = optimal_phase(m);
        plot(opt_p, metrics(m, opt_p), 'ko', 'MarkerSize', 6, 'HandleVisibility', 'off');
        if m == tx_gray + 1
            plot(opt_p, metrics(m, opt_p), 'gs', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    
    title(sprintf('Sym%d: nat=%d, gray=%d, TxBranch=%d', sym_idx, tx_nat, tx_gray, tx_gray));
    xlabel('Intra-Symbol Sample (1-16)');
    ylabel('Metric');
    grid on; xlim([1, nsps]);
    if sym_idx == 1
        legend('Location', 'northwest');
    end
end

%% ========================================================================
% 8. Summary plots: optimal phase distribution and strategy comparison
%% ========================================================================
figure('Name', 'Strategy Comparison', 'Position', [100 100 1200 400]);

% 8.1 Per-branch optimal phase distribution
subplot(1, 3, 1);
for m = 1:M
    histogram(branch_peak_pos(m, :), 0.5:1:16.5, 'FaceAlpha', 0.4, 'DisplayName', sprintf('B%d', m-1));
    hold on;
end
xline(nsps/2, 'r--', 'LineWidth', 2, 'Label', 'Nominal Mid', 'HandleVisibility', 'off');
xlabel('Peak Position (1-16)'); ylabel('Count');
title('Per-Branch Peak Position Distribution');
legend('Location', 'eastoutside'); grid on;

% 8.2 BER comparison
subplot(1, 3, 2);
bar([BER_hard_A, BER_hard_B, BER_hard_C] * 100);
set(gca, 'XTickLabel', {'Global=0', 'Global=-1', 'Per-Branch'});
ylabel('SER (%)');
title('Hard Decision SER Comparison');
grid on;

% 8.3 Discrimination ratio comparison
subplot(1, 3, 3);
bar([disc_A.mean, disc_B.mean, disc_C.mean]);
hold on;
errorbar(1:3, [disc_A.mean, disc_B.mean, disc_C.mean], ...
         [disc_A.std, disc_B.std, disc_C.std], 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Global=0', 'Global=-1', 'Per-Branch'});
ylabel('Discrimination Ratio');
title('Tx/Adjacent Ratio Comparison');
grid on;

%% ========================================================================
% 9. Conclusion output
%% ========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('1. Per-Branch Optimal Phases (relative to nominal window):\n');
for m = 1:M
    fprintf('   Branch %d (gray=%d, tone=%+d): opt_phase = %d (offset %+d from mid)\n', ...
        m-1, m-1, freq_no(m), optimal_phase(m), optimal_phase(m) - nsps/2);
end

fprintf('\n2. Performance Gain:\n');
fprintf('   Global delta=0  SER: %.4f%%\n', BER_hard_A*100);
fprintf('   Global delta=-1 SER: %.4f%%\n', BER_hard_B*100);
fprintf('   Per-Branch Opt  SER: %.4f%%\n', BER_hard_C*100);
gain = (BER_hard_A - BER_hard_C) / BER_hard_A * 100;
fprintf('   Per-Branch gain over Global=0: %.1f%% relative reduction\n', gain);

fprintf('\n3. Discrimination Ratio:\n');
fprintf('   Global delta=0:  %.3f (+-%.3f)\n', disc_A.mean, disc_A.std);
fprintf('   Global delta=-1: %.3f (+-%.3f)\n', disc_B.mean, disc_B.std);
fprintf('   Per-Branch:      %.3f (+-%.3f)\n', disc_C.mean, disc_C.std);

fprintf('\nAnalysis complete.\n');
end

%% ========================================================================
% Helper function: compute discrimination ratio
%% ========================================================================
function disc = compute_disc_ratio(bm_matrix, sym_tx, gray_enc, gry2nat)
    Nsym = length(sym_tx);
    ratios = zeros(Nsym, 1);
    for t = 1:Nsym
        tx_gray = gray_enc(sym_tx(t)+1);
        tx_val = bm_matrix(tx_gray+1, t);
        others = bm_matrix(:, t);
        others(tx_gray+1) = -inf;
        max_other = max(others);
        if max_other > 1e-6
            ratios(t) = tx_val / max_other;
        else
            ratios(t) = inf;
        end
    end
    disc.mean = mean(ratios);
    disc.std = std(ratios);
    disc.min = min(ratios);
end