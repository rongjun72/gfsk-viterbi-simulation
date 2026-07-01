function plot_viterbi_trellis()
% plot_viterbi_trellis.m
% 4-state Viterbi decoder trellis diagram visualization
% Display: 4 state nodes, 16 fully-connected transitions, forward recursion path metrics, optimal path traceback

close all; clc;

% Parameters
M = 4;              % 4 states (Gray encoding: 0,1,2,3)
T = 8;              % Show 8 symbol instants
t_states = 0:T;   % Time axis (including t=0 initial instant)

% State labels (Gray encoding)
state_labels = {'S_0 (Gray=0)', 'S_1 (Gray=1)', 'S_2 (Gray=2)', 'S_3 (Gray=3)'};
state_y = 3:-1:0;  % Y-axis: S0=3, S1=2, S2=1, S3=0 (top to bottom)

% Simulated branch metrics (example data for visualization)
% Simulated forward recursion path metrics (pm[state, t])
% Normalized pm values (randomly generated but with trend to make optimal path clear)
rng(42);
pm = zeros(M, T+1);
pm(:, 1) = [0, -2, -3, -1];  % t=1 Initial path metric (prior for prev=0)

% Forward recursion: simulate cumulative metric at each step
% Also record traceback (back[state, t] = optimal predecessor state)
back = zeros(M, T+1);

for t = 2:T+1
    for s = 1:M  % Current state (1..4 correspond to Gray 0..3)
        % Simulate 4 predecessor branches' cumulative metrics
        branch_metrics = zeros(M, 1);
        for prev = 1:M
            % Simulate branch metric: closer to correct path, higher value
            % Randomly generated but with trend added
            branch_noise = randn() * 0.5;
            % Let s=prev have slightly higher probability to get better metric (simulate some correlation)
            if prev == s
                branch_noise = branch_noise + 0.3;
            end
            branch_metrics(prev) = pm(prev, t-1) + branch_noise;
        end
        [pm(s, t), back(s, t)] = max(branch_metrics);
    end
    % Normalization (prevent overflow, consistent with code)
    pm(:, t) = pm(:, t) - max(pm(:, t));
end

%% ========================================================================
% Figure 1: Complete trellis (16 fully-connected branches)
%% ========================================================================
figure('Name', '4-State Viterbi Trellis Diagram', 'Position', [100 100 1400 700]);

% Plot state nodes
for s = 1:M
    for t = 1:T+1
        plot(t-1, state_y(s), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'LineWidth', 2);
        hold on;
        % Annotate path metric value
        if t > 1
            text(t-1, state_y(s) + 0.25, sprintf('%.2f', pm(s, t)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.2 0.2 0.6]);
        else
            text(t-1, state_y(s) + 0.25, sprintf('%.2f', pm(s, t)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.6 0.2 0.2]);
        end
    end
end

% Plot all transitions (16 per time)
% Use different transparency and color: higher branch metric, thicker/darker line
for t = 1:T
    for prev = 1:M
        for curr = 1:M
            % Calculate cumulative value of this branch (used to determine line style)
            branch_val = pm(prev, t);  % Simplified as predecessor state metric
            
            % Line thickness: according to path metric normalization
            line_width = 0.5 + 2 * (branch_val - min(pm(:,t))) / ...
                (max(pm(:,t)) - min(pm(:,t)) + 0.01);
            
            % Color: blue color scheme, darker means higher path metric
            color_intensity = (branch_val - min(pm(:,t))) / ...
                (max(pm(:,t)) - min(pm(:,t)) + 0.01);
            color_val = [0.3, 0.3 + 0.5*(1-color_intensity), 0.8];
            
            % Plot connection (with slight curve, avoid overlap)
            x = [t-1, t];
            y = [state_y(prev), state_y(curr)];
            
            % Add slight offset to avoid lines fully overlapping
            offset = (curr - prev) * 0.03;
            y_mid = (y(1) + y(2))/2 + offset;
            
            plot(x, y, 'Color', color_val, 'LineWidth', line_width, 'LineStyle', '-');
        end
    end
end

% Annotate state labels
for s = 1:M
    text(-0.5, state_y(s), state_labels{s}, ...
        'HorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
end

% Annotate time axis
for t = 0:T
    text(t, -0.8, sprintf('t=%d', t), 'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Symbol Time (t)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('State (Gray Code)', 'FontSize', 12, 'FontWeight', 'bold');
title('4-State Viterbi Trellis Diagram: 4-GFSK (16 transitions per step, all connections allowed)', ...
    'FontSize', 13, 'FontWeight', 'bold');

axis([-1 T+1 -1.5 4.5]);
set(gca, 'YTick', 0:3, 'YTickLabel', {'S_3', 'S_2', 'S_1', 'S_0'});
grid on;
box on;
hold off;

% Add legend annotation
annotation('textbox', [0.15, 0.02, 0.7, 0.05], 'String', ...
    {'Note: 4-GFSK has no state constraints (all prev→curr transitions are valid). ', ...
     'Line thickness/color intensity = path metric magnitude. ', ...
     'Numbers above nodes = normalized path metric \gamma(state, t).'}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3]);

%% ========================================================================
% Figure 2: Optimal path highlight + traceback process
%% ========================================================================
figure('Name', 'Viterbi Optimal Path Traceback', 'Position', [150 150 1400 700]);

% Replot nodes and background lines (faded)
for s = 1:M
    for t = 1:T+1
        plot(t-1, state_y(s), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'LineWidth', 2);
        hold on;
        text(t-1, state_y(s) + 0.25, sprintf('%.2f', pm(s, t)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
    end
end

% Plot all transitions (gray, faded)
for t = 1:T
    for prev = 1:M
        for curr = 1:M
            x = [t-1, t];
            y = [state_y(prev), state_y(curr)];
            plot(x, y, 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5, 'LineStyle', '-');
        end
    end
end

% Compute optimal path (traceback from t=T to t=1)
optimal_path = zeros(T+1, 1);
[~, optimal_path(end)] = max(pm(:, end));
for t = T+1:-1:2
    optimal_path(t-1) = back(optimal_path(t), t);
end

% Highlight optimal path
for t = 1:T
    s_from = optimal_path(t);
    s_to = optimal_path(t+1);
    x = [t-1, t];
    y = [state_y(s_from), state_y(s_to)];
    plot(x, y, 'r-', 'LineWidth', 3);
    
    % Annotate transfer
    mid_x = (x(1) + x(2))/2;
    mid_y = (y(1) + y(2))/2 + 0.15;
    text(mid_x, mid_y, sprintf('%d→%d', s_from-1, s_to-1), ...
        'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

% Highlight nodes on optimal path
for t = 1:T+1
    s = optimal_path(t);
    plot(t-1, state_y(s), 'ro', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    text(t-1, state_y(s) - 0.35, sprintf('S_{%d}', s-1), ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

% Annotate state labels
for s = 1:M
    text(-0.5, state_y(s), state_labels{s}, ...
        'HorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
end

for t = 0:T
    text(t, -0.8, sprintf('t=%d', t), 'HorizontalAlignment', 'center', 'FontSize', 10);
end

xlabel('Symbol Time (t)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('State (Gray Code)', 'FontSize', 12, 'FontWeight', 'bold');
title('Viterbi Traceback: Optimal Path (Red) from t=7 back to t=0', ...
    'FontSize', 13, 'FontWeight', 'bold');

axis([-1 T+1 -1.5 4.5]);
set(gca, 'YTick', 0:3, 'YTickLabel', {'S_3', 'S_2', 'S_1', 'S_0'});
grid on;
box on;
hold off;

% Add legend annotation
legend({'All states', 'Optimal path nodes', 'Optimal path transition'}, ...
    'Location', 'southwest', 'FontSize', 9);

annotation('textbox', [0.15, 0.02, 0.7, 0.05], 'String', ...
    {'Red path = Viterbi traceback result: at each step, select the predecessor with maximum cumulative metric.', ...
     'Optimal path defines the decoded Gray code sequence: S_{optimal(t)} for t=0..7.'}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3]);

%% ========================================================================
% Figure 3: Forward recursion process illustration (one-step detail)
%% ========================================================================
figure('Name', 'Viterbi Forward Recursion: One Step Detail', 'Position', [200 50 1200 500]);

% Only show one step t=2→3
focus_t = 3;  % Show recursion from t=2 to t=3

% Nodes
for s = 1:M
    plot(0, state_y(s), 'ko', 'MarkerSize', 14, 'MarkerFaceColor', [0.9 0.9 1], 'LineWidth', 2);
    hold on;
    plot(1, state_y(s), 'ko', 'MarkerSize', 14, 'MarkerFaceColor', 'w', 'LineWidth', 2);
    
    % Annotate path metric
    text(0, state_y(s) + 0.3, sprintf('\gamma_{%d}=%.2f', s-1, pm(s, focus_t-1)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.2 0.2 0.6]);
end

% All transitions (with branch metric annotation)
for prev = 1:M
    for curr = 1:M
        % SimulateBranch metric (cosine similarity ClassType)
        branch_val = randn() * 0.3 + 0.5;  % Random example value
        if prev == curr
            branch_val = branch_val + 0.2;  % Same-state transition slightly higher (simulation)
        end
        
        x = [0, 1];
        y = [state_y(prev), state_y(curr)];
        
        % Line style：According to cumulative value
        cum_val = pm(prev, focus_t-1) + branch_val;
        is_optimal = (prev == back(curr, focus_t));
        
        if is_optimal
            plot(x, y, 'g-', 'LineWidth', 2.5);
            mid_x = 0.5;
            mid_y = (y(1) + y(2))/2 + 0.1;
            text(mid_x, mid_y, sprintf('%.2f', branch_val), ...
                'Color', 'g', 'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
        else
            plot(x, y, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
        end
    end
end

% Annotate predecessor selection
for curr = 1:M
    best_prev = back(curr, focus_t);
    text(1, state_y(curr) - 0.3, sprintf('← from S_{%d}', best_prev-1), ...
        'Color', 'g', 'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    
    % New path metric
    new_pm = pm(curr, focus_t);
    text(1, state_y(curr) + 0.3, sprintf('\gamma''_{%d}=%.2f', curr-1, new_pm), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [0.2 0.6 0.2], 'FontWeight', 'bold');
end

% Annotate states
for s = 1:M
    text(-0.4, state_y(s), sprintf('S_{%d}', s-1), ...
        'HorizontalAlignment', 'right', 'FontSize', 11, 'FontWeight', 'bold');
    text(1.4, state_y(s), sprintf('S_{%d}', s-1), ...
        'HorizontalAlignment', 'left', 'FontSize', 11, 'FontWeight', 'bold');
end

text(-0.3, 4.2, 't=2 (Previous)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(1.3, 4.2, 't=3 (Current)', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('State', 'FontSize', 12, 'FontWeight', 'bold');
title({'Viterbi Forward Recursion: From t=2 to t=3', '(Green = selected predecessor, Gray = discarded)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

axis([-0.5 1.5 -1 5]);
set(gca, 'XTick', [0 1], 'XTickLabel', {'t=2', 't=3'});
set(gca, 'YTick', 0:3, 'YTickLabel', {'S_3', 'S_2', 'S_1', 'S_0'});
grid on;
box on;
hold off;

annotation('textbox', [0.1, 0.02, 0.8, 0.06], 'String', ...
    {'Forward recursion: \gamma_{curr}(t) = max_{prev}[ \gamma_{prev}(t-1) + \lambda(prev \rightarrow curr) ]', ...
     'Green arrows show the best predecessor for each current state. Gray arrows are discarded (survivor path pruning).'}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 9, ...
    'Color', [0.3 0.3 0.3]);

fprintf('\nViterbi Trellis Diagram generated.\n');
fprintf('Figure 1: Full trellis with 16 transitions per step\n');
fprintf('Figure 2: Optimal path traceback (red highlight)\n');
fprintf('Figure 3: Single-step forward recursion detail\n');

end
