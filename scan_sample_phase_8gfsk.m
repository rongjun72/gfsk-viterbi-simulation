function scan_sample_phase_8gfsk()
% scan_sample_phase_8gfsk.m
% 8-ary GFSK 采样相位偏移扫描：delta = -3, -2, -1, 0, +1, +2, +3 采样点
% 对比硬判决 BER、Viterbi BER、发送/邻支区分度
% 关键发现：
%   - delta < 0（提前采样）：Viterbi 正常（仅利用过去 ISI，可预测）
%   - delta = 0：硬判决最优
%   - delta > 0（延后采样）：Viterbi 失效（未来 ISI 泄露，不可预测）

%% 参数设置（与 gfsk_8ary_viterbi_isi 一致）
Rs=1e3; Fs=16e3; nsps=Fs/Rs; M=8; k=log2(M); h=1.0; BT=0.5; span=4; Nsym=10000;

gauss_filt=gaussdesign(BT,span,nsps); delay_gauss=grpdelay(gauss_filt,1,1)+0;
Fp=4.5e3; Fs_stop=5.5e3;
ch_filter=designfilt('lowpassfir','PassbandFrequency',Fp,'StopbandFrequency',Fs_stop, ...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',Fs);
ch_coeffs=ch_filter.Coefficients; delay_ch=grpdelay(ch_coeffs,1,1)+0;

tone_spacing=h*Rs; Fc_tone=0.75*tone_spacing;
tone_coeffs=fir1(24, Fc_tone/(Fs/2), 'low', chebwin(25, 80)); delay_tone=grpdelay(tone_coeffs,1,1)+0;

total_delay=round(delay_gauss+delay_ch+delay_tone); N_pre=ceil(total_delay/nsps)+5; N_post=N_pre; Nsym_total=Nsym+N_pre+N_post; Ns_total=Nsym_total*nsps;

gray_enc=[0;1;3;2;6;7;5;4]; gry2nat=zeros(8,1); for i=0:7, gry2nat(gray_enc(i+1)+1)=i; end
freq_no=[-7;-5;-3;-1;1;3;5;7]; tone_freq=freq_no*h*Rs/2;

fprintf('=== 8-GFSK Sample Phase Scan ===\n');

    function s=generate_gfsk(sym_seq)
        Nsym_in=length(sym_seq); Ns_in=Nsym_in*nsps;
        sym_gray=gray_enc(sym_seq+1); f_seq=freq_no(sym_gray+1);
        f_up=repelem(f_seq,nsps); f_smooth=filter(gauss_filt,1,f_up);
        dphi=2*pi*f_smooth*h*Rs/2/Fs; phase=cumsum(dphi); s=exp(1j*phase);
        if length(s)<Ns_in, s=[s;zeros(Ns_in-length(s),1)]; else, s=s(1:Ns_in); end
    end

    function bm=measure_tonemixer(r,idx)
        Ns_r=length(r); t=(0:Ns_r-1)'/Fs;
        bm=zeros(M,length(idx));
        for m=1:M
            y_mix=r.*exp(-1j*2*pi*tone_freq(m)*t);
            y_lpf=filter(tone_coeffs,1,y_mix);
            bm(m,:)=abs(y_lpf(idx)).';
        end
    end

    function det_gray=viterbi_decode(obs_matrix, ref_m)
        [~, T]=size(obs_matrix); pm=zeros(M,T); back=zeros(M,T);
        for s=1:M
            curr_g=s-1; prev_g=0;
            ref=squeeze(ref_m(prev_g+1,curr_g+1,:)); obs=obs_matrix(:,1);
            n_obs=norm(obs); n_ref=norm(ref);
            if n_obs>1e-6 && n_ref>1e-6, branch=(obs/n_obs)'*(ref/n_ref); else, branch=0; end
            pm(s,1)=branch;
        end
        for t=2:T
            for curr_g=0:M-1
                s_prime=curr_g+1; best_val=-inf; best_prev=1;
                for prev_g=0:M-1
                    s=prev_g+1;
                    ref=squeeze(ref_m(prev_g+1,curr_g+1,:)); obs=obs_matrix(:,t);
                    n_obs=norm(obs); n_ref=norm(ref);
                    if n_obs>1e-6 && n_ref>1e-6, branch=(obs/n_obs)'*(ref/n_ref); else, branch=0; end
                    val=pm(s,t-1)+branch;
                    if val>best_val, best_val=val; best_prev=s; end
                end
                pm(s_prime,t)=best_val; back(s_prime,t)=best_prev;
            end
            pm(:,t)=pm(:,t)-max(pm(:,t));
        end
        det_gray=zeros(T,1); [~,det_gray(end)]=max(pm(:,end));
        for t=T-1:-1:1, det_gray(t)=back(det_gray(t+1),t+1); end
        det_gray=det_gray-1;
    end

%% 生成固定测试信号
rng(42); sym_test=[zeros(N_pre,1); randi([0,M-1],Nsym,1); zeros(N_post,1)];
s_test=generate_gfsk(sym_test); r_test=filter(ch_coeffs,1,s_test);

%% 扫描 delta
deltas=-3:3; N_delta=length(deltas);
BER_hard=zeros(1,N_delta); BER_vit=zeros(1,N_delta); disc_ratio=zeros(1,N_delta);

for di=1:N_delta
    delta=deltas(di); fprintf('--- delta = %+d ---\n', delta);
    sample_idx=(N_pre+(0:Nsym-1))*nsps+nsps/2+total_delay+delta;
    if sample_idx(1)<1 || sample_idx(end)>Ns_total, fprintf('  SKIP (out of bounds)\n'); BER_hard(di)=NaN; BER_vit(di)=NaN; continue; end
    
    bm=measure_tonemixer(r_test,sample_idx);
    
    % 硬判决
    [~,det_gray_hard]=max(bm,[],1); det_gray_hard=det_gray_hard(:)-1; det_sym_hard=gry2nat(det_gray_hard+1);
    tx_symbols=sym_test(N_pre+1:N_pre+Nsym); tx_gray=gray_enc(tx_symbols+1);
    ratios=zeros(Nsym,1);
    for t=1:Nsym
        tx_val=bm(tx_gray(t)+1,t); others=bm(:,t); others(tx_gray(t)+1)=-inf; max_other=max(others);
        if max_other>1e-6, ratios(t)=tx_val/max_other; else, ratios(t)=inf; end
    end
    disc_ratio(di)=mean(ratios);
    BER_hard(di)=sum(det_sym_hard~=tx_symbols)/Nsym;
    
    % 预计算 ISI 参考模板
    N_guard=12; ref_metric=zeros(M,M,M);
    for prev_g=0:M-1
        for curr_g=0:M-1
            prev_nat=gry2nat(prev_g+1); curr_nat=gry2nat(curr_g+1);
            sym_seq=[zeros(N_guard,1); prev_nat; curr_nat; zeros(N_guard,1)];
            s=generate_gfsk(sym_seq); s_ch=filter(ch_coeffs,1,s);
            k_curr=N_guard+2; idx_curr=(k_curr-1)*nsps+nsps/2+total_delay+delta;
            ref_metric(prev_g+1,curr_g+1,:)=measure_tonemixer(s_ch,idx_curr);
        end
    end
    
    % Viterbi
    det_gray_vit=viterbi_decode(bm, ref_metric); det_sym_vit=gry2nat(det_gray_vit+1);
    BER_vit(di)=sum(det_sym_vit~=tx_symbols)/Nsym;
    
    fprintf('  Hard BER = %.4e | Vit BER = %.4e | Disc Ratio = %.3f\n', BER_hard(di), BER_vit(di), disc_ratio(di));
end

%% 结果汇总
fprintf('\n========== RESULT SUMMARY ==========\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Delta', 'Hard BER', 'Vit BER', 'Disc Ratio');
for di=1:N_delta
    fprintf('%+8d | %.4e | %.4e | %12.3f\n', deltas(di), BER_hard(di), BER_vit(di), disc_ratio(di));
end

[min_hard, idx_hard]=min(BER_hard); [min_vit, idx_vit]=min(BER_vit); [max_disc, idx_disc]=max(disc_ratio);
fprintf('\nBest Hard BER:   delta=%+d, BER=%.4e\n', deltas(idx_hard), min_hard);
fprintf('Best Viterbi:    delta=%+d, BER=%.4e\n', deltas(idx_vit), min_vit);
fprintf('Best Disc Ratio: delta=%+d, ratio=%.3f\n', deltas(idx_disc), max_disc);

idx0=find(deltas==0);
fprintf('\nCurrent (delta=0): Hard=%.4e, Vit=%.4e, Disc=%.3f\n', BER_hard(idx0), BER_vit(idx0), disc_ratio(idx0));

%% 绘图
figure('Name','Sample Phase Scan','Position',[100 100 1200 400]);
subplot(1,3,1); plot(deltas, BER_hard*100, 'bo-', 'LineWidth',1.5, 'MarkerSize',8); hold on; plot(deltas, BER_vit*100, 'r^-', 'LineWidth',1.5, 'MarkerSize',8); xline(0,'k--','Alpha',0.5); xlabel('Phase Offset (samples)'); ylabel('BER (%)'); title('BER vs Sampling Phase'); legend('Hard','Viterbi'); grid on;
subplot(1,3,2); semilogy(deltas, BER_hard, 'bo-', 'LineWidth',1.5, 'MarkerSize',8); hold on; semilogy(deltas, BER_vit, 'r^-', 'LineWidth',1.5, 'MarkerSize',8); xline(0,'k--','Alpha',0.5); xlabel('Phase Offset (samples)'); ylabel('BER (log)'); title('BER (log scale)'); grid on;
subplot(1,3,3); plot(deltas, disc_ratio, 'g-s', 'LineWidth',1.5, 'MarkerSize',8); hold on; plot(deltas(idx_disc), max_disc, 'ro', 'MarkerSize',12, 'LineWidth',2); xline(0,'k--','Alpha',0.5); xlabel('Phase Offset (samples)'); ylabel('Tx / Max-Adjacent Ratio'); title('Discrimination Ratio'); grid on;

fprintf('\nScan complete.\n');
end
