function test_tone_lpf_order()
% test_tone_lpf_order.m
% 测试：Tone 低通滤波器阶数对解调性能的影响（硬判决 vs ISI 感知 Viterbi）
% 扫描范围：10, 15, 20, 25, 30, 35, 40, 45, 50
% 固定参数：8-ary GFSK, h=1.0, BT=0.5, nsps=16, 无噪声 (EbN0 = 100 dB)
% 关键观察：tone LPF 阶数变化 → delay_tone 变化 → 采样位置变化 → 总 BER 变化

%% 固定参数
Rs=1e3; Fs=16e3; nsps=Fs/Rs; M=8; k=log2(M); h=1.0; BT=0.5; span=4; Nsym=5000;
tone_orders = 10:5:50;   tone_spacing = h*Rs;   Fc_tone = 0.75*tone_spacing;

% 固定信道滤波器
Fp=4.5e3; Fs_stop=5.5e3;
ch_filter = designfilt('lowpassfir', 'PassbandFrequency',Fp,'StopbandFrequency',Fs_stop, ...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',Fs);
ch_coeffs=ch_filter.Coefficients; delay_ch=grpdelay(ch_coeffs,1,1)+0;

% 固定高斯滤波器
gauss_filt=gaussdesign(BT,span,nsps); delay_gauss=grpdelay(gauss_filt,1,1)+0;

% Gray 映射
gray_enc=[0;1;3;2;6;7;5;4]; gry2nat=zeros(8,1); for i=0:7, gry2nat(gray_enc(i+1)+1)=i; end
freq_no=[-7;-5;-3;-1;1;3;5;7]; tone_freq=freq_no*h*Rs/2;

fprintf('=== Tone LPF Order Sweep Test ===\n');
fprintf('M=%d, h=%.1f, BT=%.1f, Nsym=%d, no noise\n\n', M, h, BT, Nsym);

    function s = generate_gfsk(sym_seq)
        Nsym_in=length(sym_seq); Ns_in=Nsym_in*nsps;
        sym_gray_in=gray_enc(sym_seq+1); f_seq=freq_no(sym_gray_in+1);
        f_up=repelem(f_seq,nsps); f_smooth=filter(gauss_filt,1,f_up);
        dphi=2*pi*f_smooth*h*Rs/2/Fs; phase=cumsum(dphi); s=exp(1j*phase);
        if length(s)<Ns_in, s=[s;zeros(Ns_in-length(s),1)]; else, s=s(1:Ns_in); end
    end

fprintf('Order | delay_tone | delay_total | First sample | Hard BER  | Vit BER   | Time\n');
fprintf('------|------------|-------------|--------------|-----------|-----------|------\n');

BER_hard=zeros(size(tone_orders)); BER_vit=zeros(size(tone_orders));
delay_tone_vec=zeros(size(tone_orders));

for idx=1:length(tone_orders)
    order=tone_orders(idx);
    tic_scan=tic;
    
    %% 重新设计 tone LPF
    tone_coeffs=fir1(order, Fc_tone/(Fs/2), 'low', chebwin(order+1, 80));
    delay_tone=grpdelay(tone_coeffs,1,1)+0; delay_tone_vec(idx)=delay_tone;
    
    %% 重新计算总延迟和 N_pre/N_post
    total_delay=round(delay_gauss+delay_ch+delay_tone);
    N_pre=ceil(total_delay/nsps)+5; N_post=N_pre;
    Nsym_total=Nsym+N_pre+N_post; Ns_total=Nsym_total*nsps;
    
    %% 生成固定测试信号
    rng(42); sym_tx=[zeros(N_pre,1);randi([0,M-1],Nsym,1);zeros(N_post,1)]; sym_valid=sym_tx(N_pre+1:N_pre+Nsym);
    s_raw=generate_gfsk(sym_tx); s_ch=filter(ch_coeffs,1,s_raw);
    
    %% 采样位置
    sample_idx=(N_pre+(0:Nsym-1))*nsps+nsps/2+total_delay;
    if sample_idx(1)<1 || sample_idx(end)>Ns_total
        warning('Order %d: sample_idx out of bounds, skipping', order);
        BER_hard(idx)=NaN; BER_vit(idx)=NaN; continue;
    end
    
    %% Tone-Mixer 分支度量
    Ns_r=length(s_ch); t=(0:Ns_r-1)'/Fs;
    branch_metric=zeros(M,Nsym);
    for m=1:M
        y_mix=s_ch.*exp(-1j*2*pi*tone_freq(m)*t);
        y_lpf=filter(tone_coeffs,1,y_mix);
        branch_metric(m,:)=abs(y_lpf(sample_idx)).';
    end
    
    %% 硬判决
    [~,det_gray_hard]=max(branch_metric,[],1); det_gray_hard=det_gray_hard(:)-1;
    det_sym_hard=gry2nat(det_gray_hard+1);
    bit_err_hard=0;
    for i=1:Nsym
        bit_err_hard=bit_err_hard+(bitget(sym_valid(i),3)~=bitget(det_sym_hard(i),3)) + ...
            (bitget(sym_valid(i),2)~=bitget(det_sym_hard(i),2)) + ...
            (bitget(sym_valid(i),1)~=bitget(det_sym_hard(i),1));
    end
    BER_hard(idx)=bit_err_hard/(Nsym*k);
    
    %% ISI 感知 Viterbi（参考模板含当前 tone LPF）
    N_guard=12; ref_metric=zeros(M,M,M);
    for prev_g=0:M-1
        for curr_g=0:M-1
            prev_nat=gry2nat(prev_g+1); curr_nat=gry2nat(curr_g+1);
            sym_seq=[zeros(N_guard,1);prev_nat;curr_nat;zeros(N_guard,1)];
            s=generate_gfsk(sym_seq); s_ch_ref=filter(ch_coeffs,1,s);
            k_curr=N_guard+2;
            idx_curr=(k_curr-1)*nsps+nsps/2+total_delay;
            Ns_ref=length(s_ch_ref); t_ref=(0:Ns_ref-1)'/Fs;
            bm=zeros(M,1);
            for m=1:M
                y_mix=s_ch_ref.*exp(-1j*2*pi*tone_freq(m)*t_ref);
                y_lpf=filter(tone_coeffs,1,y_mix);
                bm(m)=abs(y_lpf(idx_curr));
            end
            ref_metric(prev_g+1,curr_g+1,:)=bm;
        end
    end
    
    %% Viterbi 解码
    det_gray_vit=viterbi_decode_isi(branch_metric,ref_metric,M);
    det_sym_vit=gry2nat(det_gray_vit+1);
    bit_err_vit=0;
    for i=1:Nsym
        bit_err_vit=bit_err_vit+(bitget(sym_valid(i),3)~=bitget(det_sym_vit(i),3)) + ...
            (bitget(sym_valid(i),2)~=bitget(det_sym_vit(i),2)) + ...
            (bitget(sym_valid(i),1)~=bitget(det_sym_vit(i),1));
    end
    BER_vit(idx)=bit_err_vit/(Nsym*k);
    
    fprintf('%4d  | %8.1f   | %9d   | %12d | %.4e | %.4e | %.2f s\n', ...
        order, delay_tone, total_delay, sample_idx(1), BER_hard(idx), BER_vit(idx), toc(tic_scan));
end

%% 绘图
figure('Name','Tone LPF Order vs BER','Position',[100 100 900 500]);
semilogy(tone_orders, BER_hard, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','Hard Decision');
hold on;
semilogy(tone_orders, BER_vit, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','ISI-Aware Viterbi');
grid on; xlabel('Tone LPF Filter Order'); ylabel('Bit Error Rate (no noise)');
legend('Location','best');
title(sprintf('Tone LPF Order Sweep: 8-ary GFSK, h=%.1f, BT=%.1f, Fc=%.0f Hz, random seq, no noise', h, BT, Fc_tone));

figure('Name','Tone LPF Order vs BER + Delay','Position',[150 150 1000 500]);
yyaxis left
plot(tone_orders, BER_hard*100, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','Hard BER (%)');
hold on; plot(tone_orders, BER_vit*100, 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName','Vit BER (%)');
ylabel('Bit Error Rate (%)');
yyaxis right
plot(tone_orders, delay_tone_vec, 'g^--', 'LineWidth', 1, 'MarkerSize', 8, 'DisplayName','Tone Delay');
ylabel('Tone Filter Group Delay (samples)');
grid on; xlabel('Tone LPF Filter Order'); legend('Location','best');
title(sprintf('Tone LPF Order: BER (%%) and Group Delay | 8-ary GFSK, h=%.1f, BT=%.1f', h, BT));

fprintf('\nOptimal order (min total BER): %d\n', tone_orders(BER_vit==min(BER_vit)));
end

%% 嵌套 Viterbi 解码函数
function det_gray = viterbi_decode_isi(obs_matrix, ref_metric, M)
    [~, T] = size(obs_matrix);   pm = zeros(M, T);   back = zeros(M, T);
    for s=1:M
        curr_g=s-1; prev_g=0;
        ref=squeeze(ref_metric(prev_g+1,curr_g+1,:)); obs=obs_matrix(:,1);
        n_obs=norm(obs); n_ref=norm(ref);
        if n_obs>1e-6 && n_ref>1e-6, branch=(obs/n_obs)'*(ref/n_ref); else, branch=0; end
        pm(s,1)=branch;
    end
    for t=2:T
        for s=1:M
            best_val=-inf; best_prev=1;
            for prev=1:M
                prev_g=prev-1; curr_g=s-1;
                ref=squeeze(ref_metric(prev_g+1,curr_g+1,:)); obs=obs_matrix(:,t);
                n_obs=norm(obs); n_ref=norm(ref);
                if n_obs>1e-6 && n_ref>1e-6, branch=(obs/n_obs)'*(ref/n_ref); else, branch=0; end
                val=pm(prev,t-1)+branch;
                if val>best_val, best_val=val; best_prev=prev; end
            end
            pm(s,t)=best_val; back(s,t)=best_prev;
        end
        pm(:,t)=pm(:,t)-max(pm(:,t));
    end
    det_gray=zeros(T,1); [~,det_gray(end)]=max(pm(:,end));
    for t=T-1:-1:1, det_gray(t)=back(det_gray(t+1),t+1); end
    det_gray=det_gray-1;
end
