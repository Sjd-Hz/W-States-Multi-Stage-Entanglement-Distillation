clc; clear; close all;

% Define three ranges of r
ranges = {0.01:0.01:0.49, 0.5:0.01:0.79, 0.8:0.01:0.96};
N_max = 10000;  % Max odd N to search

for k = 1:3
    r_values = ranges{k};
    optN = zeros(size(r_values));
    optP = zeros(size(r_values));
    Ep   = zeros(size(r_values)); % efficiency = F / bestN
    
    for jj = 1:length(r_values)
        r = r_values(jj);
        F = 1 - r;
        
        bestP = 0;
        bestN = 1;
        
        % Search over odd N values
        for N = 1:2:N_max
            P1 = ((N^2 - 4*N + 3)/8) * 6 * (2*F^2/9)^2; % success prob
            P1 = min(P1,1); % cap at 1 (probabilities cannot exceed 1)
            
            if P1 >= 1
                bestP = 1;
                bestN = N;
                break; % smallest N that achieves certainty
            elseif P1 > bestP
                bestP = P1;
                bestN = N;
            end
        end
        
        optN(jj) = bestN;
        optP(jj) = bestP;
        Ep(jj)   = F / bestN; % efficiency measure
    end
    
    % Plot results for this r-range
    figure;
    subplot(2,1,1);
    plot(r_values, optN, 'ko-', 'LineWidth', 1.5);
    xlabel('Decay rate $(r)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Optimal Number of copies (N) ');
    %title([' (Range ' num2str(k) ')']);
    grid on;
    axis tight;
    
    subplot(2,1,2);
    plot(r_values, Ep, 'r-', 'LineWidth', 1.5);
   xlabel('Decay rate $(r)$', 'Interpreter', 'latex','FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold'); 
    ylabel('Optimal Efficiency');
    %title(['Optimized efficiency vs decay rate r (Range ' num2str(k) ')']);
    grid on;
    axis tight;
end
