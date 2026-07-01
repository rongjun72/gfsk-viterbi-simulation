function analyze_8gfsk_tone_metrics_intrasync()
% analyze_8gfsk_tone_metrics_intrasync.m
% 8-ary GFSK noiseless Tone-Mixer branch metric intra-symbol time-domain analysis
%
% Analysis contents:
%   1. Distribution curves of each branch metric over 16 intra-symbol samples
%   2. Whether peak positions appear at symbol midpoint (nsps/2)
%   3. Whether peak positions differ across branches
%   4. Comparison of peak positions between transmitted and interfering branches
%   5. Effect of delay on peak positions

%% ========================================================================
% 0. Parameter configuration (consistent with main simulation)
%% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 8;            % 8-ary
h       = 1.0;          % Modulation index
BT      = 0.5;          % Gaussian filter BT
span    = 4;            % Gaussian filter span

% Filter design (fully consistent with gfsk_8ary_viterbi_isi)
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

% Gray encoding
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7, gry2nat(gray_enc(i+1)+1) = i; end
freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== 8-GFSK Tone-Mixer Intra-Symbol Metric Analysis ===\n');
fprintf('Parameters: Rs=%d, Fs=%d, nsps=%d, h=%.2f, BT=%.2f, span=%d\n', Rs, Fs, nsps, h, BT, span);
fprintf('Tone spacing=%.0f Hz, Tone LPF fc=%.0f Hz (order=24)\n', tone_spacing, Fc_tone);
fprintf('Delays: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n\n', ...
    delay_gauss, delay_ch, delay_tone, total_delay);

%% ========================================================================
% 1. Helper function: generate GFSK signal
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
% 2. Generate fixed test sequence (to observe inter-symbol transitions)
%% ========================================================================
% Sequence design: covers various transition patterns
% 0,1,2,3,4,5,6,7: full symbol traversal
% 7,0,7,0: extreme frequency switching
% 0,7,1,6,2,5: large jumps
% 0,0,0,1,1,1: consecutive identical symbols
Nsym = 24;  % Short sequence for easy observation
N_pre = ceil(total_delay/nsps) + 5;
N_post = N_pre;

sym_seq = [zeros(N_pre, 1); ...
           [0; 1; 2; 3; 4; 5; 6; 7; ...   % Full traverse
            7; 0; 7; 0; ...                % Extreme switching
            0; 7; 1; 6; 2; 5; ...          % Large jump
            3; 3; 3; 4; 4; 4];             % ContinuousIdentical
           zeros(N_post, 1)];

Nsym_total = length(sym_seq);
Ns_total = Nsym_total * nsps;

fprintf('Test sequence (%d symbols): ', Nsym);
for i = N_pre+1:N_pre+Nsym
    fprintf('%d', sym_seq(i));
end
fprintf('\n');

%% ========================================================================
% 3. Generate noiseless signal and pass through channel filter
%% ========================================================================
s = generate_gfsk(sym_seq);
s_ch = filter(ch_coeffs, 1, s);

% Verify signal power
sig_pow = mean(abs(s).^2);
fprintf('Signal power = %.4f (target=1.0)\n\n', sig_pow);

%% ========================================================================
% 4. Full-time Tone-Mixer branch metrics (no decimation, keep all samples)
%% ========================================================================
t_all = (0:length(s_ch)-1)' / Fs;
branch_metric_full = zeros(M, length(s_ch));

for m = 1:M
    y_mix = s_ch .* exp(-1j * 2*pi * tone_freq(m) * t_all);
    y_lpf = filter(tone_coeffs, 1, y_mix);
    branch_metric_full(m, :) = abs(y_lpf).';  % 1 x Ns_total
end

%% ========================================================================
% 5. Analysis 1: continuous-time waveform - observe overall profile
%% ========================================================================
figure('Name', 'Branch Metrics Full Timeline', 'Position', [50 50 1400 400]);

% Select valid symbol region (skip preamble)
obs_start = (N_pre) * nsps + 1;
obs_end = (N_pre + Nsym) * nsps;

% Plot superposition of all 8 branches
hold on;
colors = lines(M);
for m = 1:M
    plot(obs_start:obs_end, branch_metric_full(m, obs_start:obs_end), ...
         'Color', colors(m,:), 'LineWidth', 1.2, 'DisplayName', sprintf('B%d', m-1));
end

% Mark symbol boundaries
for k = N_pre+1:N_pre+Nsym
    sym_start = (k-1)*nsps + 1;
    xline(sym_start, 'k--', 'Alpha', 0.3, 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
% Mark symbol midpoints
for k = N_pre+1:N_pre+Nsym
    sym_mid = (k-1)*nsps + nsps/2 + total_delay;
    xline(sym_mid, 'r:', 'Alpha', 0.5, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel('Sample Index');
ylabel('Branch Metric');
title(sprintf('8-Branch Tone-Mixer Metrics (no noise, delay=%d samples)', total_delay));
legend('Location', 'eastoutside');
grid on;
xlim([obs_start, obs_end]);

% Annotate transmitted symbols at bottom
ax = gca;
ylim_vals = ax.YLim;
for k = N_pre+1:N_pre+Nsym
    sym_start = (k-1)*nsps + 1;
    sym_mid = (k-1)*nsps + nsps/2;
    text(sym_mid, ylim_vals(1) + 0.05*(ylim_vals(2)-ylim_vals(1)), ...
         sprintf('%d', sym_seq(k)), 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end

%% ========================================================================
% 6. Analysis 2: detailed intra-symbol analysis of 16 samples per symbol
%% ========================================================================
figure('Name', 'Intra-Symbol Branch Metric Distribution', 'Position', [100 100 1400 900]);

valid_symbols = N_pre+1 : N_pre+Nsym;
n_show = min(12, Nsym);  % Show first 12 symbols

for sym_idx = 1:n_show
    k = valid_symbols(sym_idx);
    
    % This symbol's 16 sample range (delay compensation window, consistent with main simulation sampling point)
    % Symbol k at receiver side valid window = nominal start + total_delay
    start_s = (k-1)*nsps + total_delay + 1;
    end_s = (k-1)*nsps + total_delay + nsps;
    sym_samples = start_s:end_s;
    
    % ExtractThis symbol's8x16MetricMatrix
    metrics = branch_metric_full(:, sym_samples);  % 8 x 16
    
    % WhenPrevious transmitted symbol'sGrayCode
    tx_gray = gray_enc(sym_seq(k)+1);
    tx_nat = sym_seq(k);
    
    % PlotSub-Figure
    subplot(3, 4, sym_idx);
    hold on;
    
    for m = 1:M
        branch_gray = m - 1;
        if branch_gray == tx_gray
            % Transmitted branch - Red thick line
            plot(1:nsps, metrics(m,:), 'LineWidth', 2.5, 'Color', 'r', ...
                 'DisplayName', sprintf('Tx(B%d)', branch_gray));
        elseif abs(branch_gray - tx_gray) == 1 || abs(branch_gray - tx_gray) == 7
            % AdjacenttoneBranch - Orange medium line
            plot(1:nsps, metrics(m,:), 'LineWidth', 1.2, 'Color', [1 0.6 0], ...
                 'DisplayName', sprintf('Adj(B%d)', branch_gray));
        else
            % OtherBranch - Gray thin line
            plot(1:nsps, metrics(m,:), 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], ...
                 'DisplayName', sprintf('Other(B%d)', branch_gray));
        end
    end
    
    % AnnotateSymbolMidpoint（nsps/2 = 8）
    xline(nsps/2, 'b--', 'LineWidth', 1.5, 'Alpha', 0.7, 'Label', 'mid');
    
    % AnnotatePeak of transmitted branch
    [tx_peak_val, tx_peak_idx] = max(metrics(tx_gray+1, :));
    plot(tx_peak_idx, tx_peak_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
         'HandleVisibility', 'off');
    
    title(sprintf('Sym%d: nat=%d, gray=%d, TxBranch=%d', ...
                  sym_idx, tx_nat, tx_gray, tx_gray));
    xlabel('Intra-Symbol Sample (1-16)');
    ylabel('Metric');
    grid on;
    xlim([1, nsps]);
    
    % FigureExample（Only in first subplotFigureDisplay）
    if sym_idx == 1
        legend('Tx Branch', 'Adjacent', 'Other', 'Location', 'northwest');
    end
end

%% ========================================================================
% 7. Analysis 3: peak position statistics
%% ========================================================================
% Find peak position for each branch per symbol
fprintf('\n=== Peak Position Analysis (within each symbol''s 16 samples) ===\n');

all_peak_pos = zeros(M, Nsym);    % Peak value location of each branch at each symbol(1~16)
all_peak_val = zeros(M, Nsym);    % Corresponding peak value magnitude

for sym_idx = 1:Nsym
    k = valid_symbols(sym_idx);
    % Delay compensation window: symbol k at receiver side valid region
    start_s = (k-1)*nsps + total_delay + 1;
    end_s = (k-1)*nsps + total_delay + nsps;
    sym_samples = start_s:end_s;
    metrics = branch_metric_full(:, sym_samples);
    
    for m = 1:M
        [peak_val, peak_idx] = max(metrics(m,:));
        all_peak_pos(m, sym_idx) = peak_idx;
        all_peak_val(m, sym_idx) = peak_val;
    end
end

% Statistics table
fprintf('\n%-8s | %-15s | %-10s | %-12s | %-15s\n', ...
    'Branch', 'Mean Peak Pos', 'Std Pos', 'Mean Peak', 'Peak Pos Range');
fprintf('%s\n', repmat('-', 1, 75));

for m = 1:M
    mean_pos = mean(all_peak_pos(m,:));
    std_pos = std(all_peak_pos(m,:));
    mean_val = mean(all_peak_val(m,:));
    min_pos = min(all_peak_pos(m,:));
    max_pos = max(all_peak_pos(m,:));
    fprintf('%8d | %15.2f | %10.2f | %12.4f | %2d ~ %2d\n', ...
        m-1, mean_pos, std_pos, mean_val, min_pos, max_pos);
end

% Transmitted branch peak position statistics
tx_peak_pos = zeros(Nsym, 1);
for sym_idx = 1:Nsym
    k = valid_symbols(sym_idx);
    tx_gray = gray_enc(sym_seq(k)+1);
    tx_peak_pos(sym_idx) = all_peak_pos(tx_gray+1, sym_idx);
end

fprintf('\n=== Transmitted Branch Peak Position ===\n');
fprintf('Mean peak position = %.2f (target mid = %.1f)\n', mean(tx_peak_pos), nsps/2);
fprintf('Std  peak position = %.2f\n', std(tx_peak_pos));
fprintf('Range: %d ~ %d\n', min(tx_peak_pos), max(tx_peak_pos));

% Adjacent branch peak position statistics
adj_peak_pos = [];
for sym_idx = 1:Nsym
    k = valid_symbols(sym_idx);
    tx_gray = gray_enc(sym_seq(k)+1);
    for m = 1:M
        branch_gray = m - 1;
        if branch_gray ~= tx_gray && (abs(branch_gray - tx_gray) == 1 || abs(branch_gray - tx_gray) == 7)
            adj_peak_pos = [adj_peak_pos; all_peak_pos(m, sym_idx)]; %#ok<AGROW>
        end
    end
end

fprintf('\n=== Adjacent Branch Peak Position ===\n');
fprintf('Mean peak position = %.2f\n', mean(adj_peak_pos));
fprintf('Std  peak position = %.2f\n', std(adj_peak_pos));

%% ========================================================================
% 8. Analysis 4: peak position distribution histogram
%% ========================================================================
figure('Name', 'Peak Position Distribution', 'Position', [200 200 1000 400]);

% Histogram of all branches' peak positions
subplot(1, 2, 1);
all_pos = all_peak_pos(:);
histogram(all_pos, 0.5:1:16.5, 'Normalization', 'probability', 'FaceColor', [0.3 0.5 0.8]);
hold on;
xline(nsps/2, 'r--', 'LineWidth', 2, 'Label', 'Symbol Mid');
xlabel('Peak Position (within 16 samples)');
ylabel('Probability');
title('All Branches Peak Position Distribution');
grid on;
xlim([0.5, 16.5]);

% Comparison of transmitted vs adjacent branch peak positions
subplot(1, 2, 2);
histogram(tx_peak_pos, 0.5:1:16.5, 'Normalization', 'probability', ...
          'FaceColor', 'r', 'FaceAlpha', 0.6, 'DisplayName', 'Tx Branch');
hold on;
histogram(adj_peak_pos, 0.5:1:16.5, 'Normalization', 'probability', ...
          'FaceColor', [1 0.6 0], 'FaceAlpha', 0.6, 'DisplayName', 'Adjacent');
xline(nsps/2, 'b--', 'LineWidth', 2, 'Label', 'Symbol Mid');
xlabel('Peak Position (within 16 samples)');
ylabel('Probability');
title('Tx Branch vs Adjacent Branch Peak Position');
legend('Location', 'northwest');
grid on;
xlim([0.5, 16.5]);

%% ========================================================================
% 9. Analysis 5: peak positions with delay compensation
%% ========================================================================
% Actual optimal sample point = symbol midpoint + total_delay
% Within the symbol, this corresponds to a sample offset of nsps/2 + total_delay
% But since we clip at symbol boundaries, this offset may exceed the current symbol

fprintf('\n=== Delay-Compensated Analysis ===\n');
fprintf('Total delay = %d samples = %.2f symbols\n', total_delay, total_delay/nsps);

% Redefine "symbol window": from symbol start to start+16 samples, but delay-compensated
% Actual analysis: for symbol k, optimal sample = (k-1)*nsps + nsps/2 + total_delay
% This point may fall within symbol k or symbol k+1

% Re-analyze using delay-compensated "effective window"
compensated_peaks = zeros(Nsym, M);
for sym_idx = 1:Nsym
    k = valid_symbols(sym_idx);
    % Delay compensationsampling point after
    opt_sample = (k-1)*nsps + nsps/2 + total_delay;
    % Take before and after this point8sampling points（Total16）
    win_start = max(1, opt_sample - 7);
    win_end = min(Ns_total, opt_sample + 8);
    win = win_start:win_end;
    
    metrics_win = branch_metric_full(:, win);  % 8 x length(win)
    for m = 1:M
        [~, local_peak] = max(metrics_win(m,:));
        compensated_peaks(sym_idx, m) = win(local_peak);  % Global sampling index
    end
end

% Convert to offset relative to symbol midpoint
relative_offset = compensated_peaks - ((valid_symbols'-1)*nsps + nsps/2 + total_delay);

fprintf('Mean relative offset (delay-compensated window): %.2f samples\n', mean(relative_offset(:)));
fprintf('Std  relative offset: %.2f samples\n', std(relative_offset(:)));

%% ========================================================================
% 10. Analysis 6: time-varying characteristics of intra-symbol metric ratio
%% ========================================================================
% For each symbol, compute ratio of transmitted branch to second-best (as "discrimination" metric)
figure('Name', 'Tx-to-Adjacent Ratio Intra-Symbol', 'Position', [300 300 1000 400]);

ratio_matrix = zeros(Nsym, nsps);  % Discrimination degree of each symbol's 16 samples
for sym_idx = 1:Nsym
    k = valid_symbols(sym_idx);
    % Delay compensationWindow
    start_s = (k-1)*nsps + total_delay + 1;
    end_s = (k-1)*nsps + total_delay + nsps;
    sym_samples = start_s:end_s;
    metrics = branch_metric_full(:, sym_samples);  % 8x16
    
    tx_gray = gray_enc(sym_seq(k)+1);
    tx_metric = metrics(tx_gray+1, :);  % 1x16
    
    % FindEachSecond largest branch at sampling point
    second_metric = zeros(1, nsps);
    for s = 1:nsps
        sorted = sort(metrics(:, s), 'descend');
        second_metric(s) = sorted(2);
    end
    
    ratio_matrix(sym_idx, :) = tx_metric ./ second_metric;
end

% Plot discrimination curves for all symbols
hold on;
for sym_idx = 1:Nsym
    if sym_idx <= 8
        plot(1:nsps, ratio_matrix(sym_idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('Sym%d', sym_idx));
    else
        plot(1:nsps, ratio_matrix(sym_idx, :), 'LineWidth', 0.8, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end
end
xline(nsps/2, 'r--', 'LineWidth', 2, 'Label', 'Symbol Mid');
xlabel('Intra-Symbol Sample (1-16)');
ylabel('Tx / Second-Best Ratio');
title('Symbol Discrimination Ratio vs Intra-Symbol Position');
legend('Location', 'best');
grid on;
set(gca, 'YScale', 'log');

% Statistics: midpoint vs edge discrimination
mid_ratio = ratio_matrix(:, nsps/2);
edge_ratio = [ratio_matrix(:, 1); ratio_matrix(:, end)];
fprintf('\n=== Discrimination Ratio ===\n');
fprintf('At symbol mid (sample 8): mean=%.2f, std=%.2f, min=%.2f\n', ...
    mean(mid_ratio), std(mid_ratio), min(mid_ratio));
fprintf('At symbol edge (sample 1 or 16): mean=%.2f, std=%.2f, min=%.2f\n', ...
    mean(edge_ratio), std(edge_ratio), min(edge_ratio));

%% ========================================================================
% 11. Summary output
%% ========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('1. Peak Position Analysis:\n');
fprintf('   - Symbol mid point = sample %.1f\n', nsps/2);
fprintf('   - Tx branch mean peak = %.2f (offset %.2f from mid)\n', ...
    mean(tx_peak_pos), mean(tx_peak_pos) - nsps/2);
fprintf('   - Adjacent branch mean peak = %.2f (offset %.2f from mid)\n', ...
    mean(adj_peak_pos), mean(adj_peak_pos) - nsps/2);

fprintf('\n2. Key Observations:\n');
if abs(mean(tx_peak_pos) - nsps/2) < 0.5
    fprintf('   - Tx branch peaks NEAR symbol mid (good for sampling)\n');
else
    fprintf('   - Tx branch peaks OFFSET from symbol mid by %.1f samples\n', ...
        mean(tx_peak_pos) - nsps/2);
end

if std(tx_peak_pos) < 1.0
    fprintf('   - Tx branch peak position STABLE (std=%.2f)\n', std(tx_peak_pos));
else
    fprintf('   - Tx branch peak position VARIES (std=%.2f) with symbol transition\n', std(tx_peak_pos));
end

fprintf('\n3. Discrimination:\n');
fprintf('   - Best discrimination at symbol mid: ratio=%.2f\n', mean(mid_ratio));
fprintf('   - Worst discrimination at symbol edge: ratio=%.2f\n', mean(edge_ratio));
if mean(mid_ratio) > mean(edge_ratio)
    fprintf('   - Confirms: mid-point sampling is optimal\n');
end

fprintf('\nAnalysis complete.\n');
end