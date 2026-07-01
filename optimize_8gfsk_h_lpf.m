function optimize_8gfsk_h_lpf()
% optimize_8gfsk_h_lpf.m
% Under fixed BT=0.5 constraint, scan h and tone LPF parameter combinations
% Evaluate discrimination ratio and hard-decision SER under noiseless conditions
%
% Fixed parameters: Rs=1k, Fs=16k, M=8, BT=0.5, span=4
% Scan parameters:
%   h: modulation index, affects tone spacing = h*Rs and total bandwidth
%   tone LPF order: 20,24,28,32,36,40,44 (even, to avoid half-integer delay issues)
%   tone LPF fc_factor: 0.5,0.6,0.75,0.9,1.0 × tone_spacing
% Evaluation metrics:
%   1. Noiseless hard-decision SER (1000 symbols)
%   2. Tx / second-largest branch discrimination ratio (mean, min)
%   3. Total delay (samples)
%
% Constraints:
%   - Outermost tone must be within channel filter passband
%   - tone LPF order must be even (to avoid grpdelay half-integer issues)

%% ========================================================================
% 0. Fixed parameters
%% ========================================================================
Rs      = 1e3;          Fs      = 16e3;         nsps    = Fs/Rs;
M       = 8;            k       = log2(M);        BT      = 0.5;
span    = 4;            Nsym    = 1000;           Nsim    = 1;

% Gray encoding/decoding
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7, gry2nat(gray_enc(i+1)+1) = i; end
freq_no_base = [-7; -5; -3; -1; 1; 3; 5; 7];

% Scan parameters
h_values    = [0.9, 1.0, 1.1, 1.2];           % Modulation index
lpf_orders  = [20, 24, 28, 32, 36, 40, 44];   % Tone LPF order (even number)
fc_factors  = [0.5, 0.6, 0.75, 0.9, 1.0];      % fc = factor × tone_spacing

% Channel filter: fixed, passband must cover signal for all h values
% At h=1.2, outermost tone = 7×1.2×500 = 4200 Hz, plus GFSK spread ≈ ±(4200+250) = ±4450 Hz
% So Fp=5.0k, Fs_stop=6.0k
Fp = 5.0e3;   Fs_stop = 6.0e3;

fprintf('=== 8-ary GFSK Parameter Optimization (BT=0.5 fixed) ===\n');
fprintf('Fixed: Rs=%d, Fs=%d, M=%d, BT=%.1f, span=%d\n', Rs, Fs, M, BT, span);
fprintf('Scan: h=[%s], LPF_order=[%s], fc_factor=[%s]\n', ...
    num2str(h_values), num2str(lpf_orders), num2str(fc_factors));
fprintf('Total combinations: %d\n\n', length(h_values)*length(lpf_orders)*length(fc_factors));

% Result matrices
N_h = length(h_values); N_order = length(lpf_orders); N_fc = length(fc_factors);
SER_mat      = zeros(N_h, N_order, N_fc);      % Noiseless hard decision SER
disc_mean    = zeros(N_h, N_order, N_fc);      % Discrimination ratio mean value
disc_min     = zeros(N_h, N_order, N_fc);      % Discrimination ratio minimum value
delay_total  = zeros(N_h, N_order, N_fc);      % Total delay (sample count)
valid_flag   = true(N_h, N_order, N_fc);       % Is valid (tone in passband)

%% ========================================================================
% 1. Generate test signal (once, shared by all parameters)
%% ========================================================================
rng(42);
sym_test = randi([0, M-1], Nsym, 1);

%% ========================================================================
% 2. Evaluate each parameter combination
%% ========================================================================
combo_count = 0;
for hi = 1:N_h
    h = h_values(hi);
    tone_spacing = h * Rs;
    tone_freq = freq_no_base * h * Rs / 2;  % Actual frequency (Hz)
    
    % Check outermost tone whether within channel passband (include margin)
    max_tone_freq = max(abs(tone_freq));
    if max_tone_freq > Fp * 0.95
        fprintf('[h=%.2f] SKIP: max_tone=%.0f Hz > Fp=%.0f Hz (95%% margin)\n', ...
            h, max_tone_freq, Fp);
        valid_flag(hi, :, :) = false;
        continue;
    end
    
    for oi = 1:N_order
        order = lpf_orders(oi);
        
        for fi = 1:N_fc
            fc_factor = fc_factors(fi);
            Fc_tone = fc_factor * tone_spacing;
            
            combo_count = combo_count + 1;
            fprintf('[%3d/%3d] h=%.2f, order=%2d, fc=%.0f Hz (factor=%.2f) ... ', ...
                combo_count, N_h*N_order*N_fc, h, order, Fc_tone, fc_factor);
            
            % --- Filter design ---
            gauss_filt = gaussdesign(BT, span, nsps);
            delay_gauss = grpdelay(gauss_filt,1,1)+0;
            
            ch_filter = designfilt('lowpassfir', ...
                'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
                'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
                'SampleRate', Fs);
            ch_coeffs = ch_filter.Coefficients;
            delay_ch = grpdelay(ch_coeffs,1,1)+0;
            
            tone_coeffs = fir1(order, Fc_tone/(Fs/2), 'low', chebwin(order+1, 80));
            delay_tone = grpdelay(tone_coeffs,1,1)+0;
            
            total_delay = round(delay_gauss + delay_ch + delay_tone);
            N_pre  = ceil(total_delay/nsps) + 5;
            N_post = N_pre;
            Nsym_total = Nsym + N_pre + N_post;
            Ns_total = Nsym_total * nsps;
            
            % --- Generate signal ---
            sym_seq = [zeros(N_pre,1); sym_test; zeros(N_post,1)];
            
            Nsym_in = length(sym_seq);  Ns_in = Nsym_in * nsps;
            sym_gray = gray_enc(sym_seq + 1);
            f_seq = freq_no_base(sym_gray + 1);
            f_up = repelem(f_seq, nsps);
            f_smooth = filter(gauss_filt, 1, f_up);
            dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
            phase = cumsum(dphi);
            s = exp(1j * phase);
            if length(s) < Ns_in, s = [s; zeros(Ns_in - length(s), 1)];
            else, s = s(1:Ns_in); end
            
            % ChannelFiltering
            r = filter(ch_coeffs, 1, s);
            
            % --- Sampling point ---
            sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay;
            if sample_idx(1) < 1 || sample_idx(end) > Ns_total
                fprintf('SKIP (sample OOB)\n');
                valid_flag(hi, oi, fi) = false;
                continue;
            end
            
            % --- Tone-Mixer branch metrics ---
            t_all = (0:length(r)-1)' / Fs;
            branch_metric = zeros(M, Nsym);
            for m = 1:M
                y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t_all);
                y_lpf = filter(tone_coeffs, 1, y_mix);
                branch_metric(m, :) = abs(y_lpf(sample_idx)).';
            end
            
            % --- Hard decision SER ---
            [~, det_gray] = max(branch_metric, [], 1);
            det_gray = det_gray(:) - 1;
            det_sym = gry2nat(det_gray + 1);
            sym_tx = sym_test;
            SER = sum(det_sym ~= sym_tx) / Nsym;
            
            % --- Discrimination ratio ---
            ratios = zeros(Nsym, 1);
            for t = 1:Nsym
                tx_gray = gray_enc(sym_tx(t)+1);
                tx_val = branch_metric(tx_gray+1, t);
                others = branch_metric(:, t);
                others(tx_gray+1) = -inf;
                max_other = max(others);
                if max_other > 1e-6
                    ratios(t) = tx_val / max_other;
                else
                    ratios(t) = inf;
                end
            end
            
            SER_mat(hi, oi, fi)     = SER;
            disc_mean(hi, oi, fi)   = mean(ratios);
            disc_min(hi, oi, fi)    = min(ratios);
            delay_total(hi, oi, fi) = total_delay;
            
            fprintf('SER=%.3f%%, disc=%.3f (min=%.3f), delay=%d\n', ...
                SER*100, disc_mean(hi, oi, fi), disc_min(hi, oi, fi), total_delay);
        end
    end
end

%% ========================================================================
% 3. Results summary and visualization
%% ========================================================================

% 3.1 Find optimal combination
valid_idx = find(valid_flag);
[best_SER, best_idx] = min(SER_mat(valid_idx));
[best_h_i, best_o_i, best_f_i] = ind2sub(size(SER_mat), valid_idx(best_idx));

[best_disc, best_d_idx] = max(disc_mean(valid_idx));
[best_d_h, best_d_o, best_d_f] = ind2sub(size(disc_mean), valid_idx(best_d_idx));

fprintf('\n========== BEST RESULTS ==========\n');
fprintf('Best SER: h=%.2f, order=%d, fc_factor=%.2f → SER=%.4f%% (disc=%.3f, delay=%d)\n', ...
    h_values(best_h_i), lpf_orders(best_o_i), fc_factors(best_f_i), ...
    SER_mat(best_h_i, best_o_i, best_f_i)*100, ...
    disc_mean(best_h_i, best_o_i, best_f_i), ...
    delay_total(best_h_i, best_o_i, best_f_i));
fprintf('Best Disc: h=%.2f, order=%d, fc_factor=%.2f → disc=%.3f (SER=%.4f%%, delay=%d)\n', ...
    h_values(best_d_h), lpf_orders(best_d_o), fc_factors(best_d_f), ...
    disc_mean(best_d_h, best_d_o, best_d_f), ...
    SER_mat(best_d_h, best_d_o, best_d_f)*100, ...
    delay_total(best_d_h, best_d_o, best_d_f));

% Reference (h=1.0, order=24, fc=0.75)
ref_h = find(h_values == 1.0);
ref_o = find(lpf_orders == 24);
ref_f = find(fc_factors == 0.75);
if ~isempty(ref_h) && ~isempty(ref_o) && ~isempty(ref_f)
    fprintf('\nReference (h=1.0, order=24, fc=0.75):\n');
    fprintf('  SER=%.4f%%, disc=%.3f (min=%.3f), delay=%d\n', ...
        SER_mat(ref_h, ref_o, ref_f)*100, ...
        disc_mean(ref_h, ref_o, ref_f), ...
        disc_min(ref_h, ref_o, ref_f), ...
        delay_total(ref_h, ref_o, ref_f));
end

% 3.2 Visualization: SER heatmap for h=1.0, order vs fc_factor
if ~isempty(ref_h)
    figure('Name', 'SER Heatmap (h=1.0)', 'Position', [100 100 700 500]);
    ser_h1 = squeeze(SER_mat(ref_h, :, :)) * 100;  % N_order x N_fc
    imagesc(fc_factors, lpf_orders, ser_h1);
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap(jet);
    xlabel('fc factor');
    ylabel('Tone LPF Order');
    title(sprintf('Hard Decision SER (%%) - h=%.1f', h_values(ref_h)));
    
    % AnnotateEachGrid cell value
    for oi = 1:N_order
        for fi = 1:N_fc
            if ser_h1(oi, fi) < 100
                text(fc_factors(fi), lpf_orders(oi), sprintf('%.2f', ser_h1(oi, fi)), ...
                    'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 8);
            end
        end
    end
end

% 3.3 Visualization: Discrimination ratio heatmap for h=1.0, order vs fc_factor
if ~isempty(ref_h)
    figure('Name', 'Disc Ratio Heatmap (h=1.0)', 'Position', [200 200 700 500]);
    disc_h1 = squeeze(disc_mean(ref_h, :, :));
    imagesc(fc_factors, lpf_orders, disc_h1);
    set(gca, 'YDir', 'normal');
    colorbar;
    colormap(jet);
    xlabel('fc factor');
    ylabel('Tone LPF Order');
    title(sprintf('Discrimination Ratio (mean) - h=%.1f', h_values(ref_h)));
    
    for oi = 1:N_order
        for fi = 1:N_fc
            text(fc_factors(fi), lpf_orders(oi), sprintf('%.2f', disc_h1(oi, fi)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontSize', 8);
        end
    end
end

% 3.4 Optimal parameters comparison per h
figure('Name', 'Optimal per h', 'Position', [300 300 1200 400]);

subplot(1, 3, 1);
for hi = 1:N_h
    if ~any(valid_flag(hi, :, :), 'all'), continue; end
    ser_h = squeeze(SER_mat(hi, :, :));
    [min_ser, min_idx] = min(ser_h(:));
    [o_i, f_i] = ind2sub(size(ser_h), min_idx);
    bar(hi, min_ser*100, 'DisplayName', sprintf('h=%.1f, order=%d, fc=%.2f', ...
        h_values(hi), lpf_orders(o_i), fc_factors(f_i)));
    hold on;
end
set(gca, 'XTick', 1:N_h, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), h_values, 'UniformOutput', false));
ylabel('Min SER (%)');
title('Best SER per h');
legend('Location', 'best');
grid on;

subplot(1, 3, 2);
for hi = 1:N_h
    if ~any(valid_flag(hi, :, :), 'all'), continue; end
    disc_h = squeeze(disc_mean(hi, :, :));
    [max_disc, max_idx] = max(disc_h(:));
    [o_i, f_i] = ind2sub(size(disc_h), max_idx);
    bar(hi, max_disc, 'DisplayName', sprintf('h=%.1f, order=%d, fc=%.2f', ...
        h_values(hi), lpf_orders(o_i), fc_factors(f_i)));
    hold on;
end
set(gca, 'XTick', 1:N_h, 'XTickLabel', arrayfun(@(x) sprintf('%.1f', x), h_values, 'UniformOutput', false));
ylabel('Max Disc Ratio');
title('Best Disc Ratio per h');
legend('Location', 'best');
grid on;

subplot(1, 3, 3);
for hi = 1:N_h
    if ~any(valid_flag(hi, :, :), 'all'), continue; end
    ser_h = squeeze(SER_mat(hi, :, :));
    disc_h = squeeze(disc_mean(hi, :, :));
    delay_h = squeeze(delay_total(hi, :, :));
    
    % Pareto: SER vs disc ratio，Color=Delay
    valid_mask = ser_h(:) < 100;  % OnlyDisplayValidof
    scatter(disc_h(valid_mask), ser_h(valid_mask)*100, 50, delay_h(valid_mask), 'filled');
    hold on;
end
colorbar;
set(gca, 'ColorScale', 'log');
colormap(jet);
xlabel('Discrimination Ratio');
ylabel('SER (%)');
title('SER vs Disc (color=delay)');
grid on;
set(gca, 'YScale', 'log');

% 3.5 Detailed table output
fprintf('\n========== FULL RESULT TABLE (h=1.0) ==========\n');
fprintf('%-6s | %-6s | %-8s | %-10s | %-10s | %-10s | %-8s\n', ...
    'Order', 'fc_fac', 'fc(Hz)', 'SER(%)', 'Disc_mean', 'Disc_min', 'Delay');
fprintf('%s\n', repmat('-', 1, 80));
if ~isempty(ref_h)
    for oi = 1:N_order
        for fi = 1:N_fc
            fprintf('%6d | %6.2f | %8.0f | %10.4f | %10.3f | %10.3f | %8d\n', ...
                lpf_orders(oi), fc_factors(fi), fc_factors(fi)*tone_spacing, ...
                SER_mat(ref_h, oi, fi)*100, ...
                disc_mean(ref_h, oi, fi), ...
                disc_min(ref_h, oi, fi), ...
                delay_total(ref_h, oi, fi));
        end
    end
end

fprintf('\nOptimization complete.\n');
end
