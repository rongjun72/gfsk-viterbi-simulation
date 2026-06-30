function analyze_symbol_internal_metrics()
% analyze_symbol_internal_metrics.m
% 无噪声条件下：确定不同符号时 Tone-Mixer 分支度量值在符号内部 16 采样点上的分布
% 分析重点：
%   1. 各分支度量值的峰值是否位于符号中心位置
%   2. 不同分支的峰值位置是否不同
%   3. 延迟对峰值位置的影响

%% 参数设置（与主程序一致）
Rs      = 1e3;          Fs      = 16e3;         nsps    = Fs/Rs;
M       = 4;            h       = 1.0;          BT      = 0.5;          span    = 4;
Nsym_show = 8;
sym_tx = repmat([0;1;2;3], Nsym_show/4, 1);
N_guard = 10;
sym_seq = [zeros(N_guard,1); sym_tx; zeros(N_guard,1)];
Nsym_total = length(sym_seq);   Ns_total = Nsym_total * nsps;

%% 滤波器参数
gauss_filt = gaussdesign(BT, span, nsps);   delay_gauss = grpdelay(gauss_filt,1,1)+0;
Fp = 2.0e3;   Fs_stop = 2.8e3;
ch_filter = designfilt('lowpassfir', 'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, 'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;   delay_ch = grpdelay(ch_coeffs,1,1)+0;
tone_spacing = h * Rs;   Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));   delay_tone = grpdelay(tone_coeffs,1,1)+0;

total_delay = round(delay_gauss + delay_ch + delay_tone);

gray_enc = [0; 1; 3; 2];   freq_no = [-3; -1; 1; 3];   tone_freq = freq_no * h * Rs / 2;

fprintf('=== Symbol-Internal Tone-Mixer Metric Analysis ===\n');
fprintf('Total delay: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n', delay_gauss, delay_ch, delay_tone, total_delay);

%% 生成 GFSK 信号
    function s = generate_gfsk(sym_seq)
        Nsym_in = length(sym_seq);   Ns_in = Nsym_in * nsps;
        sym_gray_in = gray_enc(sym_seq + 1);
        f_seq = freq_no(sym_gray_in + 1);
        f_up = repelem(f_seq, nsps);
        f_smooth = filter(gauss_filt, 1, f_up);
        dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
        phase = cumsum(dphi);
        s = exp(1j * phase);
        if length(s) < Ns_in, s = [s; zeros(Ns_in - length(s), 1)]; else, s = s(1:Ns_in); end
    end

s = generate_gfsk(sym_seq);   s_ch = filter(ch_coeffs, 1, s);

%% 4-Branch Tone-Mixer 全时刻记录

t = (0:Ns_total-1)' / Fs;
branch_metric_full = zeros(M, Ns_total);
for m = 1:M
    y_mix = s_ch .* exp(-1j * 2*pi * tone_freq(m) * t);
    y_lpf = filter(tone_coeffs, 1, y_mix);
    branch_metric_full(m, :) = abs(y_lpf).';
end

%% 提取有效符号范围
valid_start_sample = N_guard * nsps + 1 + total_delay;
valid_end_sample = min(valid_start_sample + Nsym_show * nsps - 1, Ns_total);
metric_valid = branch_metric_full(:, valid_start_sample:valid_end_sample);
Ns_valid = size(metric_valid, 2);

%% 绘制：Figure 1 - 每个分支的符号内部分布
figure('Name', 'Symbol-Internal Branch Metrics', 'Position', [100 100 1400 900]);
for m = 1:M
    subplot(2, 2, m);   hold on;
    colors = {'b', 'r', 'g', 'm'};
    for sym_idx = 1:Nsym_show
        sym_start = (sym_idx-1) * nsps + 1;   sym_end = sym_idx * nsps;
        sym_metric = metric_valid(m, sym_start:sym_end);
        x = (0:nsps-1);
        plot(x, sym_metric, '-o', 'Color', colors{sym_idx}, 'LineWidth', 1.2, 'MarkerSize', 4, ...
             'DisplayName', sprintf('Sym%d(tx=%d)', sym_idx, sym_tx(sym_idx)));
    end
    xline(nsps/2 - 0.5, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Symbol Midpoint');
    xlabel('Sample Index within Symbol (0..15)');
    ylabel(sprintf('Branch %d Metric |y_{%d}|', m-1, m-1));
    title(sprintf('Branch %d (f = %.0f Hz)', m-1, tone_freq(m)));
    legend('Location', 'bestoutside');   grid on;   xlim([0 nsps-1]);
end
sgtitle(sprintf('Tone-Mixer Branch Metrics within Each Symbol\nh=%.1f, BT=%.1f, nsps=%d, total delay=%d', h, BT, nsps, total_delay));

%% 峰值位置统计
fprintf('\n--- Peak Position Analysis per Symbol ---\n');
peak_positions = zeros(M, Nsym_show);   peak_values = zeros(M, Nsym_show);
for sym_idx = 1:Nsym_show
    sym_start = (sym_idx-1) * nsps + 1;   sym_end = sym_idx * nsps;
    for m = 1:M
        sym_metric = metric_valid(m, sym_start:sym_end);
        [pk_val, pk_idx] = max(sym_metric);
        peak_positions(m, sym_idx) = pk_idx - 1;   peak_values(m, sym_idx) = pk_val;
    end
    fprintf('Symbol %d (tx=%d): ', sym_idx, sym_tx(sym_idx));
    for m = 1:M
        fprintf('B%d@%02d(v=%.3f) ', m-1, peak_positions(m,sym_idx), peak_values(m,sym_idx));
    end
    fprintf('\n');
end

fprintf('\n--- Summary: Peak Position Statistics ---\n');
for m = 1:M
    fprintf('Branch %d (f=%.0f Hz): mean=%.1f, std=%.1f, range=[%d,%d]\n', ...
        m-1, tone_freq(m), mean(peak_positions(m,:)), std(peak_positions(m,:)), ...
        min(peak_positions(m,:)), max(peak_positions(m,:)));
end

fprintf('\nAnalysis complete.\n');
end
