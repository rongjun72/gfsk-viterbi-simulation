function analyze_per_branch_optimal_phase_4ary()
% analyze_per_branch_optimal_phase_4ary.m
% 4-ary GFSK per-branch optimal sampling phase analysis
%
% Compare three strategies:
%   1. Global delta=0（Default）
%   2. Global delta=-1（Phase offset）
%   3. Per-Branch Optimal（EachBranchUsingEach at16Peak position within sampling point window）
%
% Differences from 8ary version:
%   - 4 tones (spacing 1000 Hz), center frequencies -1500,-500,500,1500 Hz
%   - Relaxed channel filter Fp=2.0k, Fs_stop=2.8k
%   - 36-tap tone LPF (archive default)

%% ========================================================================
% 0. Parameter configuration (consistent with gfsk_4ary_coherent_final)
%% ========================================================================
Rs      = 1e3;          Fs      = 16e3;         nsps    = Fs/Rs;
M       = 4;            k       = log2(M);        h       = 1.0;
BT      = 0.5;          span    = 4;

% Filters
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

Fp = 2.0e3;   Fs_stop = 2.8e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;   delay_ch = grpdelay(ch_coeffs,1,1)+0;

tone_spacing = h * Rs;   Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = N_pre;

% Gray encoding/decoding
gray_enc = [0; 1; 3; 2];
gry2nat = zeros(4,1);
for i = 0:3, gry2nat(gray_enc(i+1)+1) = i; end
freq_no = [-3; -1; 1; 3];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== 4-ary GFSK Per-Branch Optimal Phase Analysis ===\n');
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
% 2. Generate test signals
%% ========================================================================
rng(42);
Nsym_stat = 1000;
sym_long = [zeros(N_pre,1); randi([0,M-1],Nsym_stat,1); zeros(N_post,1)];
s_long = generate_gfsk(sym_long);
r_long = filter(ch_coeffs, 1, s_long);

sym_short = [zeros(N_pre,1); [0;1;2;3; 3;0;3;0; 0;3;1;2; 1;1;1;2;2;2]; zeros(N_post,1)];
Nsym_short = 18;
s_short = generate_gfsk(sym_short);
r_short = filter(ch_coeffs, 1, s_short);

fprintf('Generated signals: long=%d symbols, short=%d symbols\n', Nsym_stat, Nsym_short);

%% ========================================================================
% 3. Full-time branch metrics
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
% 4. Optimal phase per branch
%% ========================================================================
valid_short = N_pre+1 : N_pre+Nsym_short;
branch_peak_pos = zeros(M, Nsym_short);

for sym_idx = 1:Nsym_short
    k = valid_short(sym_idx);
    win_start = (k-1)*nsps + total_delay + 1;
    win_end = (k-1)*nsps + total_delay + nsps;
    win = win_start:win_end;
    metrics = branch_metric_full_short(:, win);
    
    for m = 1:M
        [~, peak_idx] = max(metrics(m, :));
        branch_peak_pos(m, sym_idx) = peak_idx;
    end
end

fprintf('\n--- Per-Branch Optimal Phase (relative to delay-compensated window) ---\n');
fprintf('%-8s | %-10s | %-10s | %-10s | %-10s\n', 'Branch', 'Mean Pos', 'Median', 'Mode', 'Std');
fprintf('%s\n', repmat('-', 1, 55));

optimal_phase = zeros(M, 1);
for m = 1:M
    mean_pos = mean(branch_peak_pos(m, :));
    median_pos = median(branch_peak_pos(m, :));
    mode_pos = mode(branch_peak_pos(m, :));
    std_pos = std(branch_peak_pos(m, :));
    optimal_phase(m) = round(mean_pos);
    fprintf('%8d | %10.2f | %10.1f | %10d | %10.2f\n', m-1, mean_pos, median_pos, mode_pos, std_pos);
end

%% ========================================================================
% 5. Comparison of three sampling strategies
%% ========================================================================
valid_long = N_pre+1 : N_pre+Nsym_stat;
strategies = {'Global delta=0', 'Global delta=-1', 'Per-Branch Optimal'};

sample_idx_0 = (valid_long-1)*nsps + nsps/2 + total_delay;
sample_idx_m1 = (valid_long-1)*nsps + nsps/2 + total_delay - 1;

bm_A = branch_metric_full_long(:, sample_idx_0);
bm_B = branch_metric_full_long(:, sample_idx_m1);

bm_C = zeros(M, Nsym_stat);
for m = 1:M
    sample_idx_per_branch = (valid_long-1)*nsps + total_delay + optimal_phase(m);
    bm_C(m, :) = branch_metric_full_long(m, sample_idx_per_branch);
end

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

fprintf('\n--- Hard Decision SER (no noise) ---\n');
SER_vals = [BER_hard_A, BER_hard_B, BER_hard_C];
for i = 1:3
    fprintf('%-20s: SER = %.4e (%.2f%%)\n', strategies{i}, SER_vals(i), SER_vals(i)*100);
end

%% ========================================================================
% 6. Discrimination ratio
%% ========================================================================
disc_A = compute_disc_ratio(bm_A, sym_tx, gray_enc, gry2nat);
disc_B = compute_disc_ratio(bm_B, sym_tx, gray_enc, gry2nat);
disc_C = compute_disc_ratio(bm_C, sym_tx, gray_enc, gry2nat);

fprintf('\n--- Discrimination Ratio (Tx / Max-Adjacent) ---\n');
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Global delta=0', disc_A.mean, disc_A.std, disc_A.min);
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Global delta=-1', disc_B.mean, disc_B.std, disc_B.min);
fprintf('%-20s: mean=%.3f, std=%.3f, min=%.3f\n', 'Per-Branch Optimal', disc_C.mean, disc_C.std, disc_C.min);

%% ========================================================================
% 7. Visualization: per-symbol branch curves
%% ========================================================================
figure('Name', '4-ary Per-Branch Optimal Phase', 'Position', [50 50 1400 900]);

n_show = min(8, Nsym_short);
for sym_idx = 1:n_show
    k = valid_short(sym_idx);
    win_start = (k-1)*nsps + total_delay + 1;
    win_end = (k-1)*nsps + total_delay + nsps;
    win = win_start:win_end;
    metrics = branch_metric_full_short(:, win);
    
    tx_gray = gray_enc(sym_short(k)+1);
    tx_nat = sym_short(k);
    
    subplot(4, 2, sym_idx);
    hold on;
    
    for m = 1:M
        branch_gray = m - 1;
        if branch_gray == tx_gray
            plot(1:nsps, metrics(m,:), 'LineWidth', 2.5, 'Color', 'r', 'DisplayName', sprintf('Tx(B%d)', branch_gray));
        elseif abs(branch_gray - tx_gray) == 1 || abs(branch_gray - tx_gray) == 3
            plot(1:nsps, metrics(m,:), 'LineWidth', 1.2, 'Color', [1 0.6 0], 'DisplayName', sprintf('Adj(B%d)', branch_gray));
        else
            plot(1:nsps, metrics(m,:), 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
    end
    
    xline(nsps/2, 'b--', 'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    xline(nsps/2 - 1, 'm:', 'LineWidth', 1.5, 'Alpha', 0.7, 'HandleVisibility', 'off');
    
    for m = 1:M
        opt_p = optimal_phase(m);
        plot(opt_p, metrics(m, opt_p), 'ko', 'MarkerSize', 6, 'HandleVisibility', 'off');
        if m == tx_gray + 1
            plot(opt_p, metrics(m, opt_p), 'gs', 'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end
    
    title(sprintf('Sym%d: nat=%d, gray=%d, TxBranch=%d', sym_idx, tx_nat, tx_gray, tx_gray));
    xlabel('Intra-Symbol Sample (1-16)'); ylabel('Metric');
    grid on; xlim([1, nsps]);
    if sym_idx == 1, legend('Location', 'northwest'); end
end

%% ========================================================================
% 8. Summary plots
%% ========================================================================
figure('Name', '4-ary Strategy Comparison', 'Position', [100 100 1200 400]);

subplot(1, 3, 1);
for m = 1:M
    histogram(branch_peak_pos(m, :), 0.5:1:16.5, 'FaceAlpha', 0.4, 'DisplayName', sprintf('B%d', m-1));
    hold on;
end
xline(nsps/2, 'r--', 'LineWidth', 2, 'Label', 'Nominal Mid', 'HandleVisibility', 'off');
xlabel('Peak Position (1-16)'); ylabel('Count');
title('Per-Branch Peak Position Distribution');
legend('Location', 'eastoutside'); grid on;

subplot(1, 3, 2);
bar(SER_vals * 100);
set(gca, 'XTickLabel', {'Global=0', 'Global=-1', 'Per-Branch'});
ylabel('SER (%)');
title('Hard Decision SER Comparison');
grid on;

subplot(1, 3, 3);
disc_means = [disc_A.mean, disc_B.mean, disc_C.mean];
disc_stds = [disc_A.std, disc_B.std, disc_C.std];
bar(disc_means);
hold on;
errorbar(1:3, disc_means, disc_stds, 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Global=0', 'Global=-1', 'Per-Branch'});
ylabel('Discrimination Ratio');
title('Tx/Adjacent Ratio Comparison');
grid on;

%% ========================================================================
% 9. Conclusion
%% ========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('1. Per-Branch Optimal Phases:\n');
for m = 1:M
    fprintf('   Branch %d (gray=%d, tone=%+d): opt_phase = %d (offset %+d from mid)\n', ...
        m-1, m-1, freq_no(m), optimal_phase(m), optimal_phase(m) - nsps/2);
end

fprintf('\n2. Performance:\n');
fprintf('   Global delta=0  SER: %.4f%%\n', BER_hard_A*100);
fprintf('   Global delta=-1 SER: %.4f%%\n', BER_hard_B*100);
fprintf('   Per-Branch Opt  SER: %.4f%%\n', BER_hard_C*100);

fprintf('\n3. Discrimination Ratio:\n');
fprintf('   Global delta=0:  %.3f (+-%.3f)\n', disc_A.mean, disc_A.std);
fprintf('   Global delta=-1: %.3f (+-%.3f)\n', disc_B.mean, disc_B.std);
fprintf('   Per-Branch:      %.3f (+-%.3f)\n', disc_C.mean, disc_C.std);

fprintf('\nAnalysis complete.\n');
end

%% ========================================================================
function disc = compute_disc_ratio(bm_matrix, sym_tx, gray_enc, gry2nat)
    Nsym = length(sym_tx);  M = size(bm_matrix, 1);
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
