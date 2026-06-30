function analyze_8gfsk_tone_metrics_intrasync()
% analyze_8gfsk_tone_metrics_intrasync.m
% 8-ary GFSK 无噪声条件下 Tone-Mixer 分支度量时域内部分析
% 分析发送分支与相邻分支在符号内部的峰值位置、幅度、区分度

%% 参数设置（与 gfsk_8ary_viterbi_isi.m 一致）
Rs=1e3; Fs=16e3; nsps=Fs/Rs; M=8; k=log2(M); h=1.0; BT=0.5; span=4;

% 滤波器
gauss_filt=gaussdesign(BT,span,nsps); delay_gauss=grpdelay(gauss_filt,1,1)+0;
Fp=4.5e3; Fs_stop=5.5e3;
ch_filter=designfilt('lowpassfir','PassbandFrequency',Fp,'StopbandFrequency',Fs_stop, ...
    'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',Fs);
ch_coeffs=ch_filter.Coefficients; delay_ch=grpdelay(ch_coeffs,1,1)+0;
tone_spacing=h*Rs; Fc_tone=0.75*tone_spacing;
tone_coeffs=fir1(24, Fc_tone/(Fs/2), 'low', chebwin(25, 80)); delay_tone=grpdelay(tone_coeffs,1,1)+0;

total_delay=round(delay_gauss+delay_ch+delay_tone);

gray_enc=[0;1;3;2;6;7;5;4]; gry2nat=zeros(8,1); for i=0:7, gry2nat(gray_enc(i+1)+1)=i; end
freq_no=[-7;-5;-3;-1;1;3;5;7]; tone_freq=freq_no*h*Rs/2;

fprintf('=== 8-GFSK Tone-Mixer Intra-Symbol Metric Analysis ===\n');
fprintf('Total delay=%d samples, Symbol duration=%d samples\n\n', total_delay, nsps);

    function s = generate_gfsk(sym_seq)
        Nsym_in=length(sym_seq); Ns_in=Nsym_in*nsps;
        sym_gray_in=gray_enc(sym_seq+1); f_seq=freq_no(sym_gray_in+1);
        f_up=repelem(f_seq,nsps); f_smooth=filter(gauss_filt,1,f_up);
        dphi=2*pi*f_smooth*h*Rs/2/Fs; phase=cumsum(dphi); s=exp(1j*phase);
        if length(s)<Ns_in, s=[s;zeros(Ns_in-length(s),1)]; else, s=s(1:Ns_in); end
    end

%% 测试序列：全遍历 + 跳频 + 模式切换
Nsym=24; N_pre=ceil(total_delay/nsps)+5; N_post=N_pre;
sym_seq=[zeros(N_pre,1); [0;1;2;3;4;5;6;7; 7;0;7;0; 0;7;1;6;2;5; 3;3;3;4;4;4]; zeros(N_post,1)];

s=generate_gfsk(sym_seq); s_ch=filter(ch_coeffs,1,s);

fprintf('Test sequence (%d symbols): ', Nsym);
for i=N_pre+1:N_pre+Nsym, fprintf('%d', sym_seq(i)); end, fprintf('\n');

%% 全时刻 Tone-Mixer 分支度量
Ns_r=length(s_ch); t_all=(0:Ns_r-1)'/Fs;
branch_metric_full=zeros(M, length(s_ch));
for m=1:M
    y_mix=s_ch.*exp(-1j*2*pi*tone_freq(m)*t_all);
    y_lpf=filter(tone_coeffs,1,y_mix);
    branch_metric_full(m,:)=abs(y_lpf).';
end

%% 观察窗口：有效符号范围
obs_start=(N_pre)*nsps+1; obs_end=(N_pre+Nsym)*nsps;

%% Figure 1: 全时刻波形 + 符号边界
figure('Name','Branch Metrics Full Timeline','Position',[50 50 1400 400]);
hold on; colors=lines(M);
for m=1:M
    plot(obs_start:obs_end, branch_metric_full(m,obs_start:obs_end), ...
         'Color',colors(m,:),'LineWidth',1.2,'DisplayName',sprintf('B%d',m-1));
end
for k=N_pre+1:N_pre+Nsym
    sym_start=(k-1)*nsps+1; xline(sym_start,'k--','Alpha',0.3,'LineWidth',0.8,'HandleVisibility','off');
end
for k=N_pre+1:N_pre+Nsym
    sym_mid=(k-1)*nsps+nsps/2+total_delay;
    xline(sym_mid,'r:','Alpha',0.5,'LineWidth',1.5,'HandleVisibility','off');
end
xlabel('Sample Index'); ylabel('Branch Metric');
title(sprintf('8-Branch Tone-Mixer Metrics (no noise, delay=%d samples)', total_delay));
legend('Location','eastoutside'); grid on; xlim([obs_start, obs_end]);

%% Figure 2: 每个符号内部的 16 采样点细节（前 12 个符号）
figure('Name','Intra-Symbol Branch Metric Distribution','Position',[100 100 1400 900]);
n_show=min(12,Nsym);
for sym_idx=1:n_show
    k=N_pre+sym_idx;
    win_start=(k-1)*nsps+total_delay+1; win_end=(k-1)*nsps+total_delay+nsps;
    win=win_start:win_end; metrics=branch_metric_full(:,win);
    tx_gray=gray_enc(sym_seq(k)+1); tx_nat=sym_seq(k);
    
    subplot(3,4,sym_idx); hold on;
    for m=1:M
        branch_gray=m-1;
        if branch_gray==tx_gray
            plot(1:nsps, metrics(m,:), 'LineWidth',2.5,'Color','r','DisplayName',sprintf('Tx(B%d)',branch_gray));
        elseif abs(branch_gray-tx_gray)==1 || abs(branch_gray-tx_gray)==7
            plot(1:nsps, metrics(m,:), 'LineWidth',1.2,'Color',[1 0.6 0],'DisplayName',sprintf('Adj(B%d)',branch_gray));
        else
            plot(1:nsps, metrics(m,:), 'LineWidth',0.8,'Color',[0.7 0.7 0.7],'HandleVisibility','off');
        end
    end
    xline(nsps/2,'b--','LineWidth',1.5,'Alpha',0.7,'HandleVisibility','off');
    [tx_peak_val, tx_peak_idx]=max(metrics(tx_gray+1,:));
    plot(tx_peak_idx, tx_peak_val, 'ro', 'MarkerSize',10,'LineWidth',2,'HandleVisibility','off');
    title(sprintf('Sym%d: nat=%d, gray=%d, TxBranch=%d', sym_idx, tx_nat, tx_gray, tx_gray));
    xlabel('Intra-Symbol Sample (1-16)'); ylabel('Metric'); grid on; xlim([1,nsps]);
    if sym_idx==1, legend('Location','northwest'); end
end

%% 峰值位置统计
fprintf('\n=== Peak Position Analysis (within each symbol''s 16 samples) ===\n');
all_peak_pos=zeros(M,Nsym); all_peak_val=zeros(M,Nsym);
for sym_idx=1:Nsym
    k=N_pre+sym_idx;
    start_s=(k-1)*nsps+total_delay+1; end_s=(k-1)*nsps+total_delay+nsps;
    sym_samples=start_s:end_s; metrics=branch_metric_full(:,sym_samples);
    for m=1:M
        [peak_val, peak_idx]=max(metrics(m,:));
        all_peak_pos(m,sym_idx)=peak_idx; all_peak_val(m,sym_idx)=peak_val;
    end
end

fprintf('%-8s | %-15s | %-10s | %-12s | %-15s\n', 'Branch', 'Mean Peak Pos', 'Std Pos', 'Mean Peak', 'Peak Pos Range');
fprintf('%s\n', repmat('-', 1, 75));
for m=1:M
    fprintf('%8d | %15.2f | %10.2f | %12.4f | %2d ~ %2d\n', ...
        m-1, mean(all_peak_pos(m,:)), std(all_peak_pos(m,:)), mean(all_peak_val(m,:)), ...
        min(all_peak_pos(m,:)), max(all_peak_pos(m,:)));
end

%% 发送分支峰值 vs 相邻分支峰值对比
tx_peak_pos=zeros(Nsym,1);
for sym_idx=1:Nsym
    k=N_pre+sym_idx; tx_gray=gray_enc(sym_seq(k)+1);
    tx_peak_pos(sym_idx)=all_peak_pos(tx_gray+1,sym_idx);
end

fprintf('\n=== Transmitted Branch Peak Position ===\n');
fprintf('Mean peak position = %.2f (target mid = %.1f)\n', mean(tx_peak_pos), nsps/2);
fprintf('Std  peak position = %.2f\n', std(tx_peak_pos));

adj_peak_pos=[];
for sym_idx=1:Nsym
    k=N_pre+sym_idx; tx_gray=gray_enc(sym_seq(k)+1);
    for m=1:M
        branch_gray=m-1;
        if branch_gray~=tx_gray && (abs(branch_gray-tx_gray)==1 || abs(branch_gray-tx_gray)==7)
            adj_peak_pos=[adj_peak_pos; all_peak_pos(m,sym_idx)]; %#ok<AGROW>
        end
    end
end
fprintf('Adjacent Branch Peak Position: mean=%.2f, std=%.2f\n', mean(adj_peak_pos), std(adj_peak_pos));

%% Figure 3: 峰值位置分布直方图
figure('Name','Peak Position Distribution','Position',[200 200 1000 400]);
subplot(1,2,1);
histogram(all_peak_pos(:), 0.5:1:16.5, 'Normalization','probability','FaceColor',[0.3 0.5 0.8]);
hold on; xline(nsps/2,'r--','LineWidth',2,'Label','Symbol Mid');
xlabel('Peak Position (within 16 samples)'); ylabel('Probability');
title('All Branches Peak Position Distribution'); grid on; xlim([0.5,16.5]);

subplot(1,2,2);
histogram(tx_peak_pos, 0.5:1:16.5, 'Normalization','probability','FaceColor','r','FaceAlpha',0.6,'DisplayName','Tx Branch');
hold on; histogram(adj_peak_pos, 0.5:1:16.5, 'Normalization','probability','FaceColor',[1 0.6 0],'FaceAlpha',0.6,'DisplayName','Adjacent');
xline(nsps/2,'b--','LineWidth',2,'Label','Symbol Mid','HandleVisibility','off');
xlabel('Peak Position (within 16 samples)'); ylabel('Probability');
title('Tx Branch vs Adjacent Branch Peak Position'); legend('Location','northwest'); grid on; xlim([0.5,16.5]);

fprintf('\nAnalysis complete.\n');
end
