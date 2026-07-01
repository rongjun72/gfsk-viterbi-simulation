function analyze_symbol_internal_metrics()
% analyze_symbol_internal_metrics.m
% Analyze under noiseless conditions, distribution of Tone-Mixer branch metric values over 16 intra-symbol samples for a known transmit sequence
% Core question: Do branch metric peaks appear at symbol midpoints? Do peak positions differ across branches?

%% ParametersSet（WithMain simulationConsistent）
Rs      = 1e3;          % Symbol rate (Hz)
Fs      = 16e3;         % Sampling rate (Hz)
nsps    = Fs/Rs;        % Samples per symbol = 16
M       = 4;            % 4-ary
h       = 1.0;          % Modulation index
BT      = 0.5;          % Gaussian filter BT
span    = 4;            % Gaussian filter span

% Determined transmit sequence: 0,1,2,3,0,1,2,3,... Nsym_show symbols for display
Nsym_show = 8;          % Show 8 symbols
sym_tx = repmat([0;1;2;3], Nsym_show/4, 1);  % 8 symbols

% Preamble/postamble (for filter steady-state)
N_guard = 10;
sym_seq = [zeros(N_guard,1); sym_tx; zeros(N_guard,1)];
Nsym_total = length(sym_seq);
Ns_total = Nsym_total * nsps;

%% Filter design（WithMain simulationConsistent）
gauss_filt = gaussdesign(BT, span, nsps);
delay_gauss = grpdelay(gauss_filt,1,1)+0;

% Channel filter
Fp = 2.0e3; Fs_stop = 2.8e3;
ch_filter = designfilt('lowpassfir', ...
    'PassbandFrequency', Fp, 'StopbandFrequency', Fs_stop, ...
    'PassbandRipple', 1, 'StopbandAttenuation', 80, ...
    'SampleRate', Fs);
ch_coeffs = ch_filter.Coefficients;
delay_ch = grpdelay(ch_coeffs,1,1)+0;

% Tone LPF
tone_spacing = h * Rs;
Fc_tone = 0.75 * tone_spacing;
tone_coeffs = fir1(36, Fc_tone/(Fs/2), 'low', chebwin(37, 80));
delay_tone = grpdelay(tone_coeffs,1,1)+0;

% Total delay
total_delay = round(delay_gauss + delay_ch + delay_tone);

% Gray encoding and frequency mapping (consistent with main simulation)
gray_enc = [0; 1; 3; 2];
freq_no = [-3; -1; 1; 3];
tone_freq = freq_no * h * Rs / 2;

fprintf('=== Symbol-Internal Tone-Mixer Metric Analysis ===\n');
fprintf('Total delay: gauss=%.1f, ch=%.1f, tone=%.1f, total=%d samples\n', ...
    delay_gauss, delay_ch, delay_tone, total_delay);
fprintf('Symbol duration: %d samples\n', nsps);
fprintf('Send sequence: %s\n', mat2str(sym_tx'));

%% Generate GFSK signal（Continuous-phase）
sym_gray = gray_enc(sym_seq + 1);
f_seq = freq_no(sym_gray + 1);
f_up = repelem(f_seq, nsps);
f_smooth = filter(gauss_filt, 1, f_up);
dphi = 2*pi * f_smooth * h * Rs / 2 / Fs;
phase = cumsum(dphi);
s = exp(1j * phase);

% Pass through channel filter
s_ch = filter(ch_coeffs, 1, s);

%% 4-Branch Tone-Mixer Process and recordEachSampling point'sOutput
t = (0:Ns_total-1)' / Fs;

% branch_metric_full(m, n) = magnitude of tone m at sample n
branch_metric_full = zeros(M, Ns_total);

for m = 1:M
    y_mix = s_ch .* exp(-1j * 2*pi * tone_freq(m) * t);
    y_lpf = filter(tone_coeffs, 1, y_mix);
    branch_metric_full(m, :) = abs(y_lpf).';
end

%% DetermineValidSymbolRange（ExcludePreamblePostambleofFilteringFilter transient）
% We only analyze symbols corresponding to sym_tx (remove N_guard preamble)
% Valid symbols start at index N_guard * nsps + 1
% But considering total_delay, need to shift backward
valid_start_sample = N_guard * nsps + 1 + total_delay;
valid_end_sample = valid_start_sample + Nsym_show * nsps - 1;

% Ensure no out of bounds
valid_end_sample = min(valid_end_sample, Ns_total);

% Extract full waveform of valid portion
metric_valid = branch_metric_full(:, valid_start_sample:valid_end_sample);
Ns_valid = size(metric_valid, 2);

% Time axis (in symbols, 0 at each symbol start)
t_symbol = (0:Ns_valid-1) / nsps;  % 0..Nsym_show, each symbol is 1 unit

%% AnalysisEachSymbolInternalpeak position
fprintf('\n--- Peak Position Analysis per Symbol ---\n');

% For each symbol (16 samples), find peak position for each branch (relative to symbol start)
peak_positions = zeros(M, Nsym_show);  % Branch x symbol, value as 0..15 peak value offset
peak_values = zeros(M, Nsym_show);

for sym_idx = 1:Nsym_show
    sym_start = (sym_idx-1) * nsps + 1;
    sym_end = sym_idx * nsps;
    
    for m = 1:M
        sym_metric = metric_valid(m, sym_start:sym_end);
        [pk_val, pk_idx] = max(sym_metric);
        peak_positions(m, sym_idx) = pk_idx - 1;  % 0-based: 0..15
        peak_values(m, sym_idx) = pk_val;
    end
    
    % PrintPeak of this symbolInformation
    fprintf('Symbol %d (tx=%d): ', sym_idx, sym_tx(sym_idx));
    for m = 1:M
        fprintf('B%d@%02d(v=%.3f) ', m-1, peak_positions(m,sym_idx), peak_values(m,sym_idx));
    end
    fprintf('\n');
end

fprintf('\n--- Summary: Peak Position Statistics (0=symbol start, 15=symbol end) ---\n');
for m = 1:M
    fprintf('Branch %d (f=%.0f Hz): mean=%.1f, std=%.1f, range=[%d,%d]\n', ...
        m-1, tone_freq(m), mean(peak_positions(m,:)), std(peak_positions(m,:)), ...
        min(peak_positions(m,:)), max(peak_positions(m,:)));
end

%% Plot：Figure 1 - AllSymbol superpositionDisplayEachSymbolInternalwaveform
figure('Name', 'Symbol-Internal Branch Metrics', 'Position', [100 100 1400 900]);

% Create a subplot for each branch
for m = 1:M
    subplot(2, 2, m);
    hold on;
    
    colors = {'b', 'r', 'g', 'm', 'c', 'k', [0.5 0.5 0.5], [0.8 0.4 0]};
    
    for sym_idx = 1:Nsym_show
        sym_start = (sym_idx-1) * nsps + 1;
        sym_end = sym_idx * nsps;
        sym_metric = metric_valid(m, sym_start:sym_end);
        
        % xAxis：Sampling point position within symbol 0..15
        x = (0:nsps-1);
        
        plot(x, sym_metric, '-o', 'Color', colors{sym_idx}, ...
            'LineWidth', 1.2, 'MarkerSize', 4, ...
            'DisplayName', sprintf('Sym%d(tx=%d)', sym_idx, sym_tx(sym_idx)));
    end
    
    % Mark theoretical symbolMidpoint
    xline(nsps/2 - 0.5, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Symbol Midpoint');
    
    xlabel('Sample Index within Symbol (0..15)');
    ylabel(sprintf('Branch %d Metric |y_{%d}|', m-1, m-1));
    title(sprintf('Branch %d (f = %.0f Hz)', m-1, tone_freq(m)));
    legend('Location', 'bestoutside');
    grid on;
    xlim([0 nsps-1]);
end

sgtitle(sprintf('Tone-Mixer Branch Metrics within Each Symbol (16 samples/symbol)\nh=%.1f, BT=%.1f, nsps=%d, total delay=%d', h, BT, nsps, total_delay));

%% Plot：Figure 2 - Time-continuous waveform，Mark symbolBoundaryAnd peak position
figure('Name', 'Continuous Waveform with Symbol Boundaries', 'Position', [150 150 1400 500]);

for m = 1:M
    subplot(2, 2, m);
    hold on;
    
    % PlotContinuous waveform
    t_ms = t_symbol * 1000;  % ConvertAsms
    plot(t_ms, metric_valid(m, :), 'b-', 'LineWidth', 0.8);
    
    % Mark symbolBoundary
    for sym_idx = 0:Nsym_show
        xline(sym_idx * 1000/Rs, 'g:', 'LineWidth', 1);
    end
    
    % MarkEachSymbol's peak position
    for sym_idx = 1:Nsym_show
        sym_start = (sym_idx-1) * nsps + 1;
        sym_end = sym_idx * nsps;
        sym_metric = metric_valid(m, sym_start:sym_end);
        [~, pk_idx] = max(sym_metric);
        
        % Peak position on global time axis
        pk_global = (sym_idx - 1) + (pk_idx - 1) / nsps;
        pk_val = sym_metric(pk_idx);
        
        plot(pk_global * 1000/Rs, pk_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    
    xlabel('Time (ms)');
    ylabel(sprintf('Branch %d Metric', m-1));
    title(sprintf('Branch %d (f = %.0f Hz) with Peak Markers', m-1, tone_freq(m)));
    grid on;
    xlim([0 Nsym_show * 1000/Rs]);
end

sgtitle(sprintf('Continuous Branch Metrics with Symbol Boundaries (green) and Peaks (red)\nh=%.1f, BT=%.1f, nsps=%d', h, BT, nsps));

%% Plot：Figure 3 - 3DWaterfallFigure：Symbol index vs Position within symbol vs Metric value
figure('Name', '3D Waterfall: Symbol vs Internal Position', 'Position', [200 200 1200 400]);

for m = 1:M
    subplot(1, 4, m);
    
    % BuildMatrix：Row=Symbol，Column=Position within symbol
    waterfall_data = zeros(Nsym_show, nsps);
    for sym_idx = 1:Nsym_show
        sym_start = (sym_idx-1) * nsps + 1;
        sym_end = sym_idx * nsps;
        waterfall_data(sym_idx, :) = metric_valid(m, sym_start:sym_end);
    end
    
    imagesc(0:nsps-1, 1:Nsym_show, waterfall_data);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('Sample within Symbol');
    ylabel('Symbol Index');
    title(sprintf('Branch %d (f=%.0f Hz)', m-1, tone_freq(m)));
    
    % AtEachMark peak position on symbol
    hold on;
    for sym_idx = 1:Nsym_show
        [~, pk_idx] = max(waterfall_data(sym_idx, :));
        plot(pk_idx - 1, sym_idx, 'wo', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    
    % Mark theoryMidpoint
    xline(nsps/2 - 0.5, 'w--', 'LineWidth', 2);
end

sgtitle(sprintf('Heatmap: Metric Value vs (Symbol, Position) - Red dots = Peak, White line = Midpoint\nh=%.1f, BT=%.1f, nsps=%d', h, BT, nsps));

%% Plot：Figure 4 - OnlyCorrespond toSymbol has value“Ideal”Comparison（Transmit0Only look atbranch 0，By thisClassPush）
figure('Name', 'Correct Branch vs Other Branches (per Symbol)', 'Position', [250 250 1400 500]);

% ForEachTransmitted symbol，Find whichbranchIs“Correct”of（Correspond toTransmit frequency）
correct_branch = gray_enc(sym_tx + 1) + 1;  % 1-based: which branch should be strongest

for sym_idx = 1:Nsym_show
    subplot(2, 4, sym_idx);
    hold on;
    
    sym_start = (sym_idx-1) * nsps + 1;
    sym_end = sym_idx * nsps;
    x = (0:nsps-1);
    
    % PlotAll4branches
    for m = 1:M
        sym_metric = metric_valid(m, sym_start:sym_end);
        if m == correct_branch(sym_idx)
            plot(x, sym_metric, 'r-', 'LineWidth', 2.5, 'DisplayName', sprintf('B%d (CORRECT)', m-1));
        else
            plot(x, sym_metric, 'b-', 'LineWidth', 0.8, 'DisplayName', sprintf('B%d', m-1));
        end
    end
    
    xline(nsps/2 - 0.5, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Midpoint');
    
    xlabel('Sample within Symbol');
    ylabel('Metric Value');
    title(sprintf('Sym %d (tx=%d, Gray=%d, Correct=B%d)', ...
        sym_idx, sym_tx(sym_idx), gray_enc(sym_tx(sym_idx)+1), correct_branch(sym_idx)-1));
    grid on;
    xlim([0 nsps-1]);
    
    if sym_idx == 1
        legend('Location', 'best');
    end
end

sgtitle(sprintf('Correct Branch (Red, thick) vs Others (Blue, thin) within Each Symbol\nh=%.1f, BT=%.1f, nsps=%d', h, BT, nsps));

fprintf('\nAnalysis complete.\n');

end
