function analyze_per_branch_optimal_phase_4ary()
% analyze_per_branch_optimal_phase_4ary.m
% 4-ary GFSK 每分支最优采样相位测试
% 对比：全局统一采样 vs 每分支最优相位
% 结论：每分支最优相位 INVALID，branches 必须在同一时刻比较

%% 参数设置（与 gfsk_4ary_coherent_final 一致）
Rs=1e3; Fs=16e3; nsps=Fs/Rs; M=4; k=log2(M); h=1.0; BT=0.5; span=4;

gauss_filt=gaussdesign(BT,span,nsps); delay_gauss=grpdelay(gauss_filt,1,1)+0;
Fp=2.0e3; Fs_stop=2.8e3;
ch_filter=designfilt('lowpassfir','PassbandFrequency',Fp,'StopbandFrequency',Fs_stop, ...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',Fs);
ch_coeffs=ch_filter.Coefficients; delay_ch=grpdelay(ch_coeffs,1,1)+0;

tone_spacing=h*Rs; Fc_tone=0.75*tone_spacing;
tone_coeffs=fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80)); delay_tone=grpdelay(tone_coeffs,1,1)+0;

total_delay=round(delay_gauss+delay_ch+delay_tone); N_pre=ceil(total_delay/nsps)+5; N_post=N_pre;

gray_enc=[0;1;3;2]; gry2nat=zeros(4,1); for i=0:3, gry2nat(gray_enc(i+1)+1)=i; end
freq_no=[-3;-1;1;3]; tone_freq=freq_no*h*Rs/2;

fprintf('=== Per-Branch Optimal Phase Analysis (4-ary) ===\n');

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

rng(42); Nsym_stat=1000; Nsym_short=18;
sym_long=[zeros(N_pre,1); randi([0,M-1],Nsym_stat,1); zeros(N_post,1)];
sym_short=[zeros(N_pre,1); [0;1;2;3; 3;0;3;0; 0;3;1;2; 1;1;1;2;2;2]; zeros(N_post,1)];
s_long=generate_gfsk(sym_long); r_long=filter(ch_coeffs,1,s_long);
s_short=generate_gfsk(sym_short); r_short=filter(ch_coeffs,1,s_short);

Ns_long=length(r_long); t_long=(0:Ns_long-1)'/Fs;
branch_full_long=zeros(M, Ns_long);
for m=1:M, y_mix=r_long.*exp(-1j*2*pi*tone_freq(m)*t_long); y_lpf=filter(tone_coeffs,1,y_mix); branch_full_long(m,:)=abs(y_lpf).'; end

Ns_short=length(r_short); t_short=(0:Ns_short-1)'/Fs;
branch_full_short=zeros(M, Ns_short);
for m=1:M, y_mix=r_short.*exp(-1j*2*pi*tone_freq(m)*t_short); y_lpf=filter(tone_coeffs,1,y_mix); branch_full_short(m,:)=abs(y_lpf).'; end

valid_short=N_pre+1:N_pre+Nsym_short;
branch_peak_pos=zeros(M, Nsym_short);
for sym_idx=1:Nsym_short
    k=valid_short(sym_idx);
    win_start=(k-1)*nsps+total_delay+1; win_end=(k-1)*nsps+total_delay+nsps;
    win=win_start:win_end; metrics=branch_full_short(:,win);
    for m=1:M, [~,peak_idx]=max(metrics(m,:)); branch_peak_pos(m,sym_idx)=peak_idx; end
end

optimal_phase=zeros(M,1);
for m=1:M, optimal_phase(m)=round(mean(branch_peak_pos(m,:))); end

fprintf('\n--- Per-Branch Optimal Phase ---\n');
for m=1:M, fprintf('Branch %d: mean_pos=%.2f, std=%.2f, opt=%d\n', m-1, mean(branch_peak_pos(m,:)), std(branch_peak_pos(m,:)), optimal_phase(m)); end

valid_long=N_pre+1:N_pre+Nsym_stat;
strategies={'Global delta=0', 'Global delta=-1', 'Per-Branch Optimal'};
sample_idx_0=(valid_long-1)*nsps+nsps/2+total_delay;
sample_idx_m1=(valid_long-1)*nsps+nsps/2+total_delay-1;

bm_A=branch_full_long(:,sample_idx_0);
bm_B=branch_full_long(:,sample_idx_m1);
bm_C=zeros(M, Nsym_stat);
for m=1:M
    sample_idx_per_branch=(valid_long-1)*nsps+total_delay+optimal_phase(m);
    bm_C(m,:)=branch_full_long(m, sample_idx_per_branch);
end

sym_tx=sym_long(valid_long);
BER=[0,0,0];
for t=1:Nsym_stat
    [~,dA]=max(bm_A(:,t)); dA=dA-1; sA=gry2nat(dA+1); if sA~=sym_tx(t), BER(1)=BER(1)+1; end
    [~,dB]=max(bm_B(:,t)); dB=dB-1; sB=gry2nat(dB+1); if sB~=sym_tx(t), BER(2)=BER(2)+1; end
    [~,dC]=max(bm_C(:,t)); dC=dC-1; sC=gry2nat(dC+1); if sC~=sym_tx(t), BER(3)=BER(3)+1; end
end
BER=BER/Nsym_stat;

fprintf('\n--- SER ---\n');
for i=1:3, fprintf('%-20s: SER = %.4f%%\n', strategies{i}, BER(i)*100); end

figure('Name','4-ary Per-Branch Optimal','Position',[50 50 1400 900]);
n_show=min(8,Nsym_short);
for sym_idx=1:n_show
    k=valid_short(sym_idx);
    win_start=(k-1)*nsps+total_delay+1; win_end=(k-1)*nsps+total_delay+nsps;
    win=win_start:win_end; metrics=branch_full_short(:,win);
    tx_gray=gray_enc(sym_short(k)+1); tx_nat=sym_short(k);
    
    subplot(4,2,sym_idx); hold on;
    for m=1:M
        branch_gray=m-1;
        if branch_gray==tx_gray, plot(1:nsps, metrics(m,:), 'LineWidth',2.5,'Color','r','DisplayName',sprintf('Tx(B%d)',branch_gray));
        elseif abs(branch_gray-tx_gray)==1 || abs(branch_gray-tx_gray)==3, plot(1:nsps, metrics(m,:), 'LineWidth',1.2,'Color',[1 0.6 0],'DisplayName',sprintf('Adj(B%d)',branch_gray));
        else, plot(1:nsps, metrics(m,:), 'LineWidth',0.8,'Color',[0.7 0.7 0.7],'HandleVisibility','off'); end
    end
    xline(nsps/2,'b--','LineWidth',1.5,'Alpha',0.7,'HandleVisibility','off');
    xline(nsps/2-1,'m:','LineWidth',1.5,'Alpha',0.7,'HandleVisibility','off');
    for m=1:M
        opt_p=optimal_phase(m);
        plot(opt_p, metrics(m,opt_p), 'ko', 'MarkerSize',6,'HandleVisibility','off');
        if m==tx_gray+1, plot(opt_p, metrics(m,opt_p), 'gs', 'MarkerSize',10,'LineWidth',2,'HandleVisibility','off'); end
    end
    title(sprintf('Sym%d: nat=%d, gray=%d, TxBranch=%d', sym_idx, tx_nat, tx_gray, tx_gray));
    xlabel('Intra-Symbol Sample (1-16)'); ylabel('Metric'); grid on; xlim([1,nsps]);
    if sym_idx==1, legend('Location','northwest'); end
end

figure('Name','4-ary Strategy Comparison','Position',[100 100 1200 400]);
subplot(1,3,1);
for m=1:M, histogram(branch_peak_pos(m,:), 0.5:1:16.5, 'FaceAlpha',0.4,'DisplayName',sprintf('B%d',m-1)); hold on; end
xline(nsps/2,'r--','LineWidth',2,'Label','Nominal Mid','HandleVisibility','off');
xlabel('Peak Position (1-16)'); ylabel('Count'); title('Per-Branch Peak Position'); legend('Location','eastoutside'); grid on;

subplot(1,3,2); bar(BER*100); set(gca,'XTickLabel',{'Global=0','Global=-1','Per-Branch'}); ylabel('SER (%)'); title('Hard Decision SER'); grid on;

subplot(1,3,3);
[dA_mean,dA_std,dA_min]=compute_disc(bm_A,sym_tx,gray_enc,gry2nat);
[dB_mean,dB_std,dB_min]=compute_disc(bm_B,sym_tx,gray_enc,gry2nat);
[dC_mean,dC_std,dC_min]=compute_disc(bm_C,sym_tx,gray_enc,gry2nat);
bar([dA_mean,dB_mean,dC_mean]); hold on; errorbar(1:3,[dA_mean,dB_mean,dC_mean],[dA_std,dB_std,dC_std],'k.','LineWidth',1.5);
set(gca,'XTickLabel',{'Global=0','Global=-1','Per-Branch'}); ylabel('Discrimination Ratio'); title('Tx/Adjacent Ratio'); grid on;

fprintf('\n=== CONCLUSION ===\n');
fprintf('Per-Branch Optimal Phase is INVALID for Viterbi/argmax detection.\n');
fprintf('Global delta=0  SER: %.4f%%\n', BER(1)*100);
fprintf('Global delta=-1 SER: %.4f%%\n', BER(2)*100);
fprintf('Per-Branch Opt  SER: %.4f%%\n', BER(3)*100);

end

function [mean_d, std_d, min_d] = compute_disc(bm_matrix, sym_tx, gray_enc, gry2nat)
    Nsym=length(sym_tx); ratios=zeros(Nsym,1);
    for t=1:Nsym
        tx_gray=gray_enc(sym_tx(t)+1); tx_val=bm_matrix(tx_gray+1,t);
        others=bm_matrix(:,t); others(tx_gray+1)=-inf; max_other=max(others);
        if max_other>1e-6, ratios(t)=tx_val/max_other; else, ratios(t)=inf; end
    end
    mean_d=mean(ratios); std_d=std(ratios); min_d=min(ratios);
end
