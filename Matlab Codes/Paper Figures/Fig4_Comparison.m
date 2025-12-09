clc
clear all;
close all;

% Define parameters
r_values = 0.01:0.005:0.99;
w_values = 0.01:0.005:0.99;

% Preallocate arrays for efficiency
n_r = length(r_values);
E_optimal = NaN(1, n_r);  % Use NaN for uncomputed values
W_optimal = NaN(1, n_r);
M_optimal = NaN(1, n_r);
E_mine5 = zeros(1, n_r);  % Preallocate

% Main loop over r values
for jj = 1:n_r
    r = r_values(jj);
    F = 1 - r;
    d = r;
    
    E_mine5(jj) = (24 * F^5) / 405;  % Calculate for all r
    
    if r < 0.25
        % Preallocate for inner loop
        E_values = zeros(1, length(w_values));
        m_values = zeros(1, length(w_values));
        
        for ii = 1:length(w_values)
            w = w_values(ii);
            
            % Calculate m
            l0 = (1 - w) * (d / (1 - d));
            tt = log((3 * 1e-6) / (1 - 1e-6)) / log(3 * l0);
            m_val = real(ceil(log2(tt)));
            
            if ~isfinite(m_val)
                m_val = 8;  % Default value if Inf
            end
            m_values(ii) = m_val;
            
            % Calculate E
            E_values(ii) = calculate_E_m_N(3, m_val, w, d);
        end
        
        % Find optimal values
        [max_E, max_idx] = max(E_values);
        W_optimal(jj) = w_values(max_idx);
        M_optimal(jj) = m_values(max_idx);
        E_optimal(jj) = calculate_E_m_N(3, M_optimal(jj), W_optimal(jj), d);
    end
end

% Plot results
figure(1);
plot(r_values, E_mine5, 'r-', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
hold on;
plot(r_values, E_optimal, 'b--*', 'LineWidth', 1,'MarkerSize', 4, 'MarkerFaceColor', 'b'); 
xlabel('Decay rate $(r)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Efficiency $(E)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
 %Create the legend and set its properties
leg = legend('$E_{\mathrm{ours}}^5$','$E_{\mathrm{WM}}$' , 'Location', 'best');
set(leg, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
axis tight;
hold off;

figure(2);
plot(r_values, W_optimal, 'b--*', 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Decay rate $(r)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Optimal Weak measurement strength $(q)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
% Create the legend and set its properties
%leg = legend('$q_{\mathrm{opt}}$', 'Location', 'best');
%set(leg, 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
axis tight;