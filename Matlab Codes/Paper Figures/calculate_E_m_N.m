function E_m_N = calculate_E_m_N(N, m, w, d)
    % Calculate E_m^{(N)} with d_bar = 1 - d and w_bar = 1 - w
    % E_m^{(N)} = p_w * ∏(p_i/2) from i=1 to m
    % = w_bar^{N-1} * d_bar * (1 + lambda_m) * ∏(1/(2N(1+lambda_i))) from i=0 to m-1
    
    % Inputs:
    %   N: parameter N
    %   m: parameter m
    %   w: parameter w (used to calculate w_bar = 1 - w)
    %   d: parameter d (used to calculate d_bar = 1 - d)
    
    % Calculate derived parameters
    w_bar = 1 - w;
    d_bar = 1 - d;
    
    % Calculate lambda_i for i = 0, 1, ..., m
    lambda = zeros(1, m+1);
    for i = 0:m
        lambda(i+1) = (1/N) * (N * w_bar*(d/d_bar))^(2^i);
    end
    
    % Calculate the product term: ∏(1/(2N(1+lambda_i))) from i=0 to m-1
    product_term = 1;
    for i = 0:(m-1)
        product_term = product_term * (1 / (2 * N * (1 + lambda(i+1))));
    end
    
    % Calculate E_m^{(N)}
    E_m_N = (w_bar^(N-1)) * d_bar * (1 + lambda(m+1)) * product_term;
    
    % Display results
    %fprintf('Parameters:\n');
    %fprintf('  N = %d, m = %d\n', N, m);
    %fprintf('  w = %.4f, w_bar = 1 - w = %.4f\n', w, w_bar);
    %fprintf('  d = %.4f, d_bar = 1 - d = %.4f\n', d, d_bar);
    %fprintf('lambda values for i = 0 to %d:\n', m);
    %for i = 0:m
    %    fprintf('  lambda_%d = %.6e\n', i, lambda(i+1));
    %end
    %fprintf('Product term = %.6e\n', product_term);
   % fprintf('E_%d^{(%d)} = %.6e\n\n', m, N, E_m_N);
end

