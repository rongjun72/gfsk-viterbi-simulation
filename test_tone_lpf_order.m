function test_tone_lpf_order()
% test_tone_lpf_order.m
% Test: Effect of tone mixer lowpass filter order on noiseless BER (hard decision vs ISI-aware Viterbi)
% Scan orders: 10, 15, 20, 25, 30, 35, 40, 45, 50
% Fixed parameters: 8-ary GFSK, h=1.0, BT=0.5, nsps=16, noiseless (EbN0 = 100 dB)
% Key observation: tone LPF order change → delay_tone change → sample index change → detection BER change
%               Also LPF order affects adjacent channel isolation, thus affecting hard-decision branch leakage
% Since Viterbi reference templates also depend on delay_tone, each order requires recomputing reference templates.
%
% Note: This test uses random Tx sequence + N_pre=ceil(total_delay/nsps)+5,
%     fully consistent with gfsk_8ary_viterbi_isi.m to ensure comparable results.

%% ========================================================================
% 0. Fixed parameters
% ========================================================================
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 8;            % 8-ary
k       = log2(M);      % 3 bits/symbol
h       = 1.0;          % Modulation index
BT      = 0.5;          % Gaussian filtering BT
span    = 4;            % Gaussian filtering span
Nsym    = 5000;         % Valid symbol count (reduced to accelerate scan)

% Scanned order range
tone_orders = 10:5:50;  % [10, 15, 20, 25, 30, 35, 40, 45, 50]

% Fixed tone LPF cutoff frequency (consistent with main script)
tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;  % 750 Hz

% Fixed channel filter
Fp = 4.5e3;   Fs_stop = 5.5e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;
delay_ch = grpdelay(ch_coeffs,1,1)+0;

% Fixed Gaussian filter
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

% Gray encoding (8-ary)
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7
    gry2nat(gray_enc(i+1)+1) = i;
end
freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== Tone LPF Order Sweep Test ===\n');
fprintf('M=%d, h=%.1f, BT=%.1f, Nsym=%d, no noise\n\n', M, h, BT, Nsym);

%% ========================================================================
% 1. Helper function: generate GFSK signal (noiseless)
% ========================================================================
    function s = generate_gfsk(sym_seq)
        Nsym_in = length(sym_seq);
        sym_gray_in = gray_enc(sym_seq + 1);
        f_seq = freq_no(sym_gray_in + 1);
        f_up = repelem(f_seq, nsps);
        f_smooth = filter(gauss_filt, 1, f_up);
        dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
        phase = cumsum(dphi);
        s = exp(1j * phase);
        Ns_in = Nsym_in * nsps;
        if length(s) < Ns_in
            s = [s; zeros(Ns_in - length(s), 1)];
        else
            s = s(1:Ns_in);
        end
    end

%% ========================================================================
% 2. Scan orders
%% ========================================================================
BER_hard = zeros(size(tone_orders));
BER_vit  = zeros(size(tone_orders));
delay_tone_vec = zeros(size(tone_orders));

fprintf('Order | delay_tone | delay_total | First sample | Hard BER  | Vit BER   | Time\n');
fprintf('------|------------|-------------|--------------|-----------|-----------|------\n');

for idx = 1:length(tone_orders)
    order = tone_orders(idx);
    tic_scan = tic;
    
    %% 2.1 Redesign tone LPF
    tone_coeffs = fir1(order, Fc_tone/(Fs/2), 'low', chebwin(order+1, 80));
    delay_tone = grpdelay(tone_coeffs,1,1)+0;
    delay_tone_vec(idx) = delay_tone;
    
    %% 2.2 Re-CalculateTotal delayAnd N_pre/N_post（With gfsk_8ary_viterbi_isi.m Consistent）
    total_delay = round(delay_gauss + delay_ch + delay_tone);
    N_pre = ceil(total_delay/nsps) + 5;
    N_post = ceil(total_delay/nsps) + 5;
    
    Nsym_total = Nsym + N_pre + N_post;
    Ns_total = Nsym_total * nsps;
    
    %% 2.3 Generate random transmit sequence (consistent with gfsk_8ary_viterbi_isi.m, rng(42) fixed)
    rng(42);  % Fixed seed, ensure reproducibility
    sym_tx = [zeros(N_pre, 1); randi([0, M-1], Nsym, 1); zeros(N_post, 1)];
    sym_valid = sym_tx(N_pre+1 : N_pre+Nsym);
    
    %% 2.4 Generate GFSK Signal passing through channelFiltering
    s_raw = generate_gfsk(sym_tx);
    s_ch = filter(ch_coeffs, 1, s_raw);  % Through channel filtering (noiseless)
    
    %% 2.5 CalculateSampling index
    sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay;
    
    % SafeCheck
    if sample_idx(1) < 1 || sample_idx(end) > Ns_total
        warning('Order %d: sample_idx out of bounds, skipping', order);
        BER_hard(idx) = NaN; BER_vit(idx) = NaN;
        continue;
    end
    
    %% 2.6 Measure tone-mixer Branch metric
    Ns_r = length(s_ch);
    t = (0:Ns_r-1)' / Fs;
    branch_metric = zeros(M, Nsym);
    for m = 1:M
        y_mix = s_ch .* exp(-1j * 2*pi * tone_freq(m) * t);
        y_lpf = filter(tone_coeffs, 1, y_mix);
        branch_metric(m, :) = abs(y_lpf(sample_idx)).';
    end
    
    %% 2.7 Hard decision
    [~, det_gray_hard] = max(branch_metric, [], 1);
    det_gray_hard = det_gray_hard(:) - 1;
    det_sym_hard = gry2nat(det_gray_hard + 1);
    
    % Hard decision BER（Bit-level）
    bit_err_hard = 0;
    for i = 1:Nsym
        bit_err_hard = bit_err_hard + ...
            (bitget(sym_valid(i), 3) ~= bitget(det_sym_hard(i), 3)) + ...
            (bitget(sym_valid(i), 2) ~= bitget(det_sym_hard(i), 2)) + ...
            (bitget(sym_valid(i), 1) ~= bitget(det_sym_hard(i), 1));
    end
    BER_hard(idx) = bit_err_hard / (Nsym * k);
    
    %% 2.8 Re-Calculate Viterbi Reference template（ISI Aware）
    N_guard = 12;
    ref_metric = zeros(M, M, M);  % (prev_gray, curr_gray, branch_idx)
    
    for prev_g = 0:M-1
        for curr_g = 0:M-1
            prev_nat = gry2nat(prev_g + 1);
            curr_nat = gry2nat(curr_g + 1);
            sym_seq = [zeros(N_guard, 1); prev_nat; curr_nat; zeros(N_guard, 1)];
            s = generate_gfsk(sym_seq);
            s_ch_ref = filter(ch_coeffs, 1, s);
            
            k_curr = N_guard + 2;
            % UsingWhenPrevious order's total_delay CalculateSampling point
            idx_curr = (k_curr - 1) * nsps + nsps/2 + total_delay;
            
            Ns_ref = length(s_ch_ref);
            t_ref = (0:Ns_ref-1)' / Fs;
            bm = zeros(M, 1);
            for m = 1:M
                y_mix = s_ch_ref .* exp(-1j * 2*pi * tone_freq(m) * t_ref);
                y_lpf = filter(tone_coeffs, 1, y_mix);
                bm(m) = abs(y_lpf(idx_curr));
            end
            ref_metric(prev_g+1, curr_g+1, :) = bm;
        end
    end
    
    %% 2.9 Viterbi Decode
    det_gray_vit = viterbi_decode_isi(branch_metric, ref_metric, M);
    det_sym_vit = gry2nat(det_gray_vit + 1);
    
    bit_err_vit = 0;
    for i = 1:Nsym
        bit_err_vit = bit_err_vit + ...
            (bitget(sym_valid(i), 3) ~= bitget(det_sym_vit(i), 3)) + ...
            (bitget(sym_valid(i), 2) ~= bitget(det_sym_vit(i), 2)) + ...
            (bitget(sym_valid(i), 1) ~= bitget(det_sym_vit(i), 1));
    end
    BER_vit(idx) = bit_err_vit / (Nsym * k);
    
    fprintf('%4d  | %8.1f   | %9d   | %12d | %.4e | %.4e | %.2f s\n', ...
        order, delay_tone, total_delay, sample_idx(1), BER_hard(idx), BER_vit(idx), toc(tic_scan));
end

%% ========================================================================
% 3. Visualization results
% ========================================================================

% Figure 1: Order vs noiseless BER (linear scale)
figure('Name', 'Tone LPF Order vs BER (No Noise, Random Sequence)', 'Position', [100 100 900 500]);
semilogy(tone_orders, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard Decision');
hold on;
semilogy(tone_orders, BER_vit, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'ISI-Aware Viterbi');
grid on;
xlabel('Tone LPF Filter Order (tap count - 1)');
ylabel('Bit Error Rate (no noise)');
legend('Location', 'best');
title(sprintf('Tone LPF Order Sweep: 8-ary GFSK, h=%.1f, BT=%.1f, Fc=%.0f Hz, random seq, no noise', h, BT, Fc_tone));

% Figure 2: Order vs noiseless BER + delay info
figure('Name', 'Tone LPF Order vs BER + Delay', 'Position', [150 150 1000 500]);

yyaxis left
plot(tone_orders, BER_hard*100, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Hard BER (%)');
hold on;
plot(tone_orders, BER_vit*100, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Vit BER (%)');
ylabel('Bit Error Rate (%)');
valid_ber = [BER_hard(:); BER_vit(:)];
valid_ber = valid_ber(~isnan(valid_ber));
if isempty(valid_ber) || all(valid_ber == 0)
    ymax = 0.01;  % at least 1% for plotting
else
    ymax = max(valid_ber);
end
ylim([0 ymax*100*1.2]);

yyaxis right
plot(tone_orders, delay_tone_vec, 'g^--', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName', 'Tone Delay (samples)');
ylabel('Tone Filter Group Delay (samples)');

grid on;
xlabel('Tone LPF Filter Order');
legend('Location', 'best');
title(sprintf('Tone LPF Order: BER (%%) and Group Delay | 8-ary GFSK, h=%.1f, BT=%.1f, random seq', h, BT));

% Figure 3: Viterbi gain (dB) vs order
figure('Name', 'Viterbi Gain vs Tone LPF Order', 'Position', [200 200 600 400]);
valid = BER_vit > 0 & BER_hard > 0;
gain_dB = zeros(size(tone_orders));
gain_dB(valid) = 10*log10(BER_hard(valid) ./ BER_vit(valid));
gain_dB(~valid) = 0;

plot(tone_orders, gain_dB, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Tone LPF Filter Order');
ylabel('Viterbi Gain (dB) = 10*log_{10}(BER_{hard}/BER_{vit})');
title('ISI-Aware Viterbi Gain over Hard Decision (No Noise, Random Sequence)');

%% ========================================================================
% 4. Results summary
%% ========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Order | delay_tone |  Hard BER   |  Vit BER   | Gain(dB)\n');
fprintf('------|------------|-------------|------------|----------\n');
for i = 1:length(tone_orders)
    if BER_vit(i) > 0 && BER_hard(i) > 0
        g = 10*log10(BER_hard(i)/BER_vit(i));
        fprintf('%4d  | %8.1f   | %.4e  | %.4e | %8.2f\n', ...
            tone_orders(i), delay_tone_vec(i), BER_hard(i), BER_vit(i), g);
    else
        fprintf('%4d  | %8.1f   | %.4e  | %.4e |     N/A\n', ...
            tone_orders(i), delay_tone_vec(i), BER_hard(i), BER_vit(i));
    end
end

fprintf('\nOptimal order (min total BER): %d\n', tone_orders(BER_vit == min(BER_vit)));
valid_gain = gain_dB > 0;
if any(valid_gain)
    fprintf('Optimal order (max Viterbi gain): %d\n', tone_orders(gain_dB == max(gain_dB(valid_gain))));
end

end

%% ========================================================================
% Nested function: Viterbi decoding (8-state ISI-aware)
% Input parameter ref_metric varies with order, passed as parameter
%% ========================================================================
function det_gray = viterbi_decode_isi(obs_matrix, ref_metric, M)
    [M_v, T] = size(obs_matrix);
    
    pm = zeros(M_v, T);
    back = zeros(M_v, T);
    
    % t=1: initialize using prev=0 (preamble) ISI reference
    for s = 1:M_v
        curr_g = s-1;
        prev_g = 0;
        ref = squeeze(ref_metric(prev_g+1, curr_g+1, :));
        obs = obs_matrix(:, 1);
        
        n_obs = norm(obs);
        n_ref = norm(ref);
        if n_obs > 1e-6 && n_ref > 1e-6
            obs_n = obs / n_obs;
            ref_n = ref / n_ref;
            branch = obs_n' * ref_n;  % cosine similarity
        else
            branch = 0;
        end
        
        pm(s, 1) = branch;
    end
    
    % t=2:T
    for t = 2:T
        for s = 1:M_v
            best_val = -inf;
            best_prev = 1;
            for prev = 1:M_v
                prev_g = prev-1;
                curr_g = s-1;
                ref = squeeze(ref_metric(prev_g+1, curr_g+1, :));
                obs = obs_matrix(:, t);
                
                n_obs = norm(obs);
                n_ref = norm(ref);
                if n_obs > 1e-6 && n_ref > 1e-6
                    obs_n = obs / n_obs;
                    ref_n = ref / n_ref;
                    branch = obs_n' * ref_n;
                else
                    branch = 0;
                end
                
                val = pm(prev, t-1) + branch;
                if val > best_val
                    best_val = val;
                    best_prev = prev;
                end
            end
            pm(s, t) = best_val;
            back(s, t) = best_prev;
        end
        pm(:, t) = pm(:, t) - max(pm(:, t));
    end
    
    % Traceback
    det_gray = zeros(T, 1);
    [~, det_gray(end)] = max(pm(:, end));
    for t = T-1:-1:1
        det_gray(t) = back(det_gray(t+1), t+1);
    end
    det_gray = det_gray - 1;
end
