function optimize_8gfsk_h_lpf()
% optimize_8gfsk_h_lpf.m
% BT=0.5 固定约束下，联合优化 h 和 tone LPF 参数
% 扫描参数：
%   h: [0.9, 1.0, 1.1, 1.2]
%   tone LPF order: [20,24,28,32,36,40,44]（偶数阶）
%   fc_factor: [0.5,0.6,0.75,0.9,1.0] × tone_spacing
% 评价指标：
%   1. 无噪声硬判决 SER（1000 符号）
%   2. 发送/邻支区分度比（均值、最小值）
%   3. 总延迟
% 关键结论：
%   - h=1.2, order=20, fc=1.0 → SER=3.5%（最优）
%   - h=1.0, order=20, fc=1.0 → SER=6.1%
%   - 延迟主导：低阶 + 宽 fc = 更好采样位置

%% 固定参数
Rs=1e3; Fs=16e3; nsps=Fs/Rs; M=8; k=log2(M); BT=0.5; span=4; Nsym=1000; Nsim=1;

gray_enc=[0;1;3;2;6;7;5;4]; gry2nat=zeros(8,1); for i=0:7, gry2nat(gray_enc(i+1)+1)=i; end
freq_no_base=[-7;-5;-3;-1;1;3;5;7];

h_values=[0.9, 1.0, 1.1, 1.2];
lpf_orders=[20,24,28,32,36,40,44];
fc_factors=[0.5,0.6,0.75,0.9,1.0];

% 信道滤波器（固定，通带需覆盖 h=1.2 最大 tone ≈ 4200 Hz）
Fp=5.0e3; Fs_stop=6.0e3;

fprintf('=== 8-ary GFSK Parameter Optimization (BT=0.5 fixed) ===\n');
fprintf('Fixed: Rs=%d, Fs=%d, M=%d, BT=%.1f, span=%d\n', Rs, Fs, M, BT, span);
fprintf('Scan: h=%s, orders=%s, fc_factors=%s\n', num2str(h_values), num2str(lpf_orders), num2str(fc_factors));
fprintf('Total combinations: %d\n\n', length(h_values)*length(lpf_orders)*length(fc_factors));

N_h=length(h_values); N_order=length(lpf_orders); N_fc=length(fc_factors);
SER_mat=zeros(N_h,N_order,N_fc); disc_mean=zeros(N_h,N_order,N_fc); disc_min=zeros(N_h,N_order,N_fc); delay_total=zeros(N_h,N_order,N_fc); valid_flag=true(N_h,N_order,N_fc);

rng(42); sym_test=randi([0,M-1],Nsym,1);

combo_count=0;
for hi=1:N_h
    h=h_values(hi); tone_spacing=h*Rs; tone_freq=freq_no_base*h*Rs/2;
    max_tone_freq=max(abs(tone_freq));
    if max_tone_freq>Fp*0.95, fprintf('[h=%.2f] SKIP: max_tone=%.0f Hz > Fp*0.95=%.0f Hz\n', h, max_tone_freq, Fp*0.95); valid_flag(hi,:,:)=false; continue; end
    
    for oi=1:N_order
        order=lpf_orders(oi);
        for fi=1:N_fc
            fc_factor=fc_factors(fi); Fc_tone=fc_factor*tone_spacing;
            combo_count=combo_count+1;
            fprintf('[%3d/%3d] h=%.2f, order=%2d, fc=%.0f Hz (factor=%.2f) ... ', combo_count, N_h*N_order*N_fc, h, order, Fc_tone, fc_factor);
            
            % 滤波器设计
            gauss_filt=gaussdesign(BT,span,nsps); delay_gauss=grpdelay(gauss_filt,1,1)+0;
            ch_filter=designfilt('lowpassfir','PassbandFrequency',Fp,'StopbandFrequency',Fs_stop, ...
                'PassbandRipple',1,'StopbandAttenuation',80,'SampleRate',Fs);
            ch_coeffs=ch_filter.Coefficients; delay_ch=grpdelay(ch_coeffs,1,1)+0;
            tone_coeffs=fir1(order, Fc_tone/(Fs/2), 'low', chebwin(order+1, 80)); delay_tone=grpdelay(tone_coeffs,1,1)+0;
            
            total_delay=round(delay_gauss+delay_ch+delay_tone);
            N_pre=ceil(total_delay/nsps)+5; N_post=N_pre; Nsym_total=Nsym+N_pre+N_post; Ns_total=Nsym_total*nsps;
            
            % 信号生成
            sym_seq=[zeros(N_pre,1); sym_test; zeros(N_post,1)];
            Nsym_in=length(sym_seq); Ns_in=Nsym_in*nsps;
            sym_gray=gray_enc(sym_seq+1); f_seq=freq_no_base(sym_gray+1);
            f_up=repelem(f_seq,nsps); f_smooth=filter(gauss_filt,1,f_up);
            dphi=2*pi*f_smooth*h*Rs/2/Fs; phase=cumsum(dphi); s=exp(1j*phase);
            if length(s)<Ns_in, s=[s;zeros(Ns_in-length(s),1)]; else, s=s(1:Ns_in); end
            
            r=filter(ch_coeffs,1,s);
            
            % 采样
            sample_idx=(N_pre+(0:Nsym-1))*nsps+nsps/2+total_delay;
            if sample_idx(1)<1 || sample_idx(end)>Ns_total, fprintf('SKIP (sample OOB)\n'); valid_flag(hi,oi,fi)=false; continue; end
            
            % Tone-Mixer 分支度量
            Ns_r=length(r); t_all=(0:Ns_r-1)'/Fs;
            branch_metric=zeros(M,Nsym);
            for m=1:M
                y_mix=r.*exp(-1j*2*pi*tone_freq(m)*t_all);
                y_lpf=filter(tone_coeffs,1,y_mix);
                branch_metric(m,:)=abs(y_lpf(sample_idx)).';
            end
            
            % 硬判决 SER
            [~,det_gray]=max(branch_metric,[],1); det_gray=det_gray(:)-1; det_sym=gry2nat(det_gray+1);
            SER=sum(det_sym~=sym_test)/Nsym;
            
            % 区分度比
            ratios=zeros(Nsym,1);
            for t=1:Nsym
                tx_gray=gray_enc(sym_test(t)+1); tx_val=branch_metric(tx_gray+1,t);
                others=branch_metric(:,t); others(tx_gray+1)=-inf; max_other=max(others);
                if max_other>1e-6, ratios(t)=tx_val/max_other; else, ratios(t)=inf; end
            end
            
            SER_mat(hi,oi,fi)=SER; disc_mean(hi,oi,fi)=mean(ratios); disc_min(hi,oi,fi)=min(ratios); delay_total(hi,oi,fi)=total_delay;
            fprintf('SER=%.3f%%, disc=%.3f (min=%.3f), delay=%d\n', SER*100, mean(ratios), min(ratios), total_delay);
        end
    end
end

%% 最佳结果
valid_idx=find(valid_flag);
[best_SER, best_idx]=min(SER_mat(valid_idx));
[best_h_i, best_o_i, best_f_i]=ind2sub(size(SER_mat), valid_idx(best_idx));

[best_disc, best_d_idx]=max(disc_mean(valid_idx));
[best_d_h, best_d_o, best_d_f]=ind2sub(size(disc_mean), valid_idx(best_d_idx));

fprintf('\n========== BEST RESULTS ==========\n');
fprintf('Best SER: h=%.2f, order=%d, fc_factor=%.2f → SER=%.4f%% (disc=%.3f, delay=%d)\n', ...
    h_values(best_h_i), lpf_orders(best_o_i), fc_factors(best_f_i), SER_mat(best_h_i,best_o_i,best_f_i)*100, ...
    disc_mean(best_h_i,best_o_i,best_f_i), delay_total(best_h_i,best_o_i,best_f_i));
fprintf('Best Disc: h=%.2f, order=%d, fc_factor=%.2f → disc=%.3f (SER=%.4f%%, delay=%d)\n', ...
    h_values(best_d_h), lpf_orders(best_d_o), fc_factors(best_d_f), disc_mean(best_d_h,best_d_o,best_d_f), ...
    SER_mat(best_d_h,best_d_o,best_d_f)*100, delay_total(best_d_h,best_d_o,best_d_f));

ref_h=find(h_values==1.0); ref_o=find(lpf_orders==24); ref_f=find(fc_factors==0.75);
if ~isempty(ref_h)&&~isempty(ref_o)&&~isempty(ref_f)
    fprintf('\nReference (h=1.0, order=24, fc=0.75):\n');
    fprintf('  SER=%.4f%%, disc=%.3f (min=%.3f), delay=%d\n', ...
        SER_mat(ref_h,ref_o,ref_f)*100, disc_mean(ref_h,ref_o,ref_f), disc_min(ref_h,ref_o,ref_f), delay_total(ref_h,ref_o,ref_f));
end

%% 绘图
if ~isempty(ref_h)
    figure('Name','SER Heatmap (h=1.0)','Position',[100 100 700 500]);
    ser_h1=squeeze(SER_mat(ref_h,:,:))*100;
    imagesc(fc_factors, lpf_orders, ser_h1); set(gca,'YDir','normal'); colorbar; colormap(jet);
    xlabel('fc factor'); ylabel('Tone LPF Order'); title(sprintf('Hard Decision SER (%%) - h=%.1f', h_values(ref_h)));
    for oi=1:N_order, for fi=1:N_fc, if ser_h1(oi,fi)<100, text(fc_factors(fi),lpf_orders(oi),sprintf('%.2f',ser_h1(oi,fi)),'HorizontalAlignment','center','Color','w','FontSize',8); end, end, end
    
    figure('Name','Disc Ratio Heatmap (h=1.0)','Position',[200 200 700 500]);
    disc_h1=squeeze(disc_mean(ref_h,:,:));
    imagesc(fc_factors, lpf_orders, disc_h1); set(gca,'YDir','normal'); colorbar; colormap(jet);
    xlabel('fc factor'); ylabel('Tone LPF Order'); title(sprintf('Discrimination Ratio (mean) - h=%.1f', h_values(ref_h)));
    for oi=1:N_order, for fi=1:N_fc, text(fc_factors(fi),lpf_orders(oi),sprintf('%.2f',disc_h1(oi,fi)),'HorizontalAlignment','center','Color','w','FontSize',8); end, end
end

figure('Name','Optimal per h','Position',[300 300 1200 400]);
subplot(1,3,1);
for hi=1:N_h
    if ~any(valid_flag(hi,:,:),'all'), continue; end
    ser_h=squeeze(SER_mat(hi,:,:)); [min_ser, min_idx]=min(ser_h(:)); [o_i,f_i]=ind2sub(size(ser_h),min_idx);
    bar(hi, min_ser*100, 'DisplayName', sprintf('h=%.1f, order=%d, fc=%.2f', h_values(hi), lpf_orders(o_i), fc_factors(f_i))); hold on;
end
set(gca,'XTick',1:N_h,'XTickLabel',arrayfun(@(x) sprintf('%.1f',x), h_values, 'UniformOutput', false));
ylabel('Min SER (%)'); title('Best SER per h'); legend('Location','best'); grid on;

subplot(1,3,2);
for hi=1:N_h
    if ~any(valid_flag(hi,:,:),'all'), continue; end
    disc_h=squeeze(disc_mean(hi,:,:)); [max_disc, max_idx]=max(disc_h(:)); [o_i,f_i]=ind2sub(size(disc_h),max_idx);
    bar(hi, max_disc, 'DisplayName', sprintf('h=%.1f, order=%d, fc=%.2f', h_values(hi), lpf_orders(o_i), fc_factors(f_i))); hold on;
end
set(gca,'XTick',1:N_h,'XTickLabel',arrayfun(@(x) sprintf('%.1f',x), h_values, 'UniformOutput', false));
ylabel('Max Disc Ratio'); title('Best Disc Ratio per h'); legend('Location','best'); grid on;

subplot(1,3,3);
for hi=1:N_h
    if ~any(valid_flag(hi,:,:),'all'), continue; end
    ser_h=squeeze(SER_mat(hi,:,:)); disc_h=squeeze(disc_mean(hi,:,:)); delay_h=squeeze(delay_total(hi,:,:));
    valid_mask=ser_h(:)<100; scatter(disc_h(valid_mask), ser_h(valid_mask)*100, 50, delay_h(valid_mask), 'filled'); hold on;
end
colorbar; colormap(jet); xlabel('Discrimination Ratio'); ylabel('SER (%)'); title('SER vs Disc (color=delay)'); grid on; set(gca,'YScale','log');

fprintf('\nOptimization complete.\n');
end
