function gfsk_8ary_coherent_final()
% gfsk_8ary_coherent_final.m
% 8-ary GFSK 相干解调仿真 含 仿真报告版
% 功能：tone-mixer 软检测、理论BER对比、无SNR下误码率

%% ========================================================================
% 0. 仿真参数
% ========================================================================
Rs      = 1e3;          % 符号率 (Hz)
Fs      = 16e3;         % 采样率 (Hz)
nsps    = Fs/Rs;        % 每符号采样点数 = 16
M       = 8;            % 8进制
k       = log2(M);      % 3 bits/symbol
h       = 1.0;          % 调制指数，相邻tone间距 = h*Rs = 1000 Hz
BT      = 0.5;          % 高斯滤波BT
span    = 4;            % 高斯滤波span（符号数）
Nsym    = 10000;        % 有效符号数

EbN0_dB = 12*log10(1:1.9:20)/log10(20);  % 对数EbN0分布0~12dB
Nsim    = 1;            % 每EbN0点Monte Carlo次数

RUN_FLOOR_ANALYSIS = true;
RUN_H_SCAN         = true;

%% ========================================================================
% 1. 滤波器参数与延迟计算
% ========================================================================
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

% 信道滤波器
Fp = 4.5e3;   Fs_stop = 5.5e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;
delay_ch = grpdelay(ch_filter.Coefficients,1,1)+0;

% Tone低通滤波器
 tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% 总延迟
total_delay = round(delay_gauss + delay_ch + delay_tone);
N_pre  = ceil(total_delay/nsps) + 5;
N_post = ceil(total_delay/nsps) + 5;
Nsym_total = Nsym + N_pre + N_post;
Ns_total   = Nsym_total * nsps;

sample_idx = (N_pre + (0:Nsym-1)) * nsps + nsps/2 + total_delay;
if sample_idx(1) < 1 || sample_idx(end) > Ns_total
    error('采样点越界');
end

fprintf('=== 8-ary GFSK Coherent Demodulation ===\n');

%% ========================================================================
% 2. Gray 映射表
% ========================================================================
gray_enc = [0; 1; 3; 2; 6; 7; 5; 4];
gry2nat = zeros(8,1);
for i = 0:7
    gry2nat(gray_enc(i+1)+1) = i;
end

freq_no = [-7; -5; -3; -1; 1; 3; 5; 7];
tone_freq = freq_no * h * Rs / 2;

%% ========================================================================
% 3. 信号生成函数
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
% 4. 检测函数
%% ========================================================================
    function [det_sym, det_gray, branch_metric] = detect_coherent(r, Nsym_in, sample_idx_in)
        Ns_r = length(r);
        t = (0:Ns_r-1)' / Fs;
        branch_metric = zeros(M, Nsym_in);
        for m = 1:M
            y_mix = r .* exp(-1j * 2*pi * tone_freq(m) * t);
            y_lpf = filter(tone_coeffs, 1, y_mix);
            branch_metric(m, :) = abs(y_lpf(sample_idx_in)).';
        end
        [~, det_gray] = max(branch_metric, [], 1);
        det_gray = (det_gray(:) - 1);
        det_sym = gry2nat(det_gray + 1);
    end

% ... 完整代码（约500行）包含仿真循环、理论BER计算、绘图 ...

end
