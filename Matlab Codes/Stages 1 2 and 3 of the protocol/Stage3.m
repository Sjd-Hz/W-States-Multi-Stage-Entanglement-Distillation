clc; clear;

%% Symbolic Setup
syms F real
assume(F >= 0 & F <= 1);

%% Basic Qubit States and Operators
H = [1; 0];    % |0⟩
V = [0; 1];    % |1⟩
Id2 = eye(2);  
Id4 = eye(4);
Id8 = eye(8);
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];

%% Original 3-Qubit W State
psi_W = (1/sqrt(3)) * (kron(V, kron(H, H)) + kron(H, kron(V, H)) + kron(H, kron(H, V)));
rho_in0 =F*(psi_W * psi_W')+(1-F)*((kron(H, kron(H, H)))*(kron(H, kron(H, H)) )');
pin1=1/sqrt(4)*(kron(V, kron(H, H)) + kron(H, kron(V, H)) + kron(H, kron(H, V)) + kron(V, kron(V, V))); 
rin1=pin1*pin1';
rho_W=psi_W * psi_W';

%% initial state 
rho_total = kron(rho_in0, rin1)

%% Display Non-Zero Diagonal Elements
diag_elems = diag(rho_total);
fprintf('Non-zero diagonal elements:\n');
for idx = 1:length(diag_elems)
    if ~isequal(diag_elems(idx), sym(0))
        fprintf('rho(%d,%d) = %s\n', idx, idx, char(diag_elems(idx)));
    end
end

%% Apply Local CNOTs (Purification ⟶ Sacrificial)
Rc= CNOTn(6,1,4)*CNOTn(6,2,5)*CNOTn(6,3,6)*rho_total*(CNOTn(6,1,4)*CNOTn(6,2,5)*CNOTn(6,3,6))';

% Define computational basis projectors for sacrificial qubits
proj0 = [1, 0; 0, 0];  % |0⟩⟨0|
proj1 = [0, 0; 0, 1];  % |1⟩⟨1|


% Prepare basis combinations from 000 to 110
outcomes = dec2bin(0:7, 3);  % 000 to 110

for idx = 1:size(outcomes, 1)
    outcome = outcomes(idx, :);  % e.g., '000'

    % Build projector Mz for the current measurement outcome
    Mz = Id8;  % Start with identity for system qubits
    for k = 1:3
        if outcome(k) == '0'
            Mz = kron(Mz, proj0);
        else
            Mz = kron(Mz, proj1);
        end
    end

    % Apply projection and partial trace sacrificial qubits (4,5,6)
    rho_filtered = Mz * Rc * Mz;
    rho_purified = PartialTrace(rho_filtered, [4, 5, 6]);

    % Normalize
    p_success = simplify(trace(rho_purified));
    if p_success ~= 0
        rho_norm = simplify(rho_purified / p_success);
    else
        rho_norm = zeros(size(rho_purified));
    end

    % Display results
    fprintf('\nMeasurement outcome: |%s⟩\n', outcome);
    fprintf('Success probability: %s\n', char(p_success));
    disp('Final normalized state:');
    disp(rho_norm);
    %disp(rho_purified);
end
