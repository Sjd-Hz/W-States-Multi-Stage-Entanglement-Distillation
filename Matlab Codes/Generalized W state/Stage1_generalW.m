% Protocol 1(F*W*W'+(1-F)*000)
clc; clear;

%% Symbolic Setup
syms F  f real
assume(0 <= F <= 1);
assume(0 <= f <= 1);
%% Basic Qubit States and Operators
H = [1; 0];    % |0⟩
V = [0; 1];    % |1⟩
Id2 = eye(2);  
Id4 = eye(4);
Id8 = eye(8);
X = [0 1; 1 0];
Y = [1 0; 0 -1];
X3=kron(X,kron(X,X));
%% Original 3-Qubit W State
psi_W = (1/2) * (kron(V, kron(H, H)) + kron(H, kron(V, H)) +sqrt(2)* kron(H, kron(H, V)));
rho_inF = F*(psi_W * psi_W')+(1-F)* (kron(H, kron(H, H)))* (kron(H, kron(H, H)))';
rho_inf = f*(psi_W * psi_W')+(1-f)* (kron(H, kron(H, H)))* (kron(H, kron(H, H)))';
rw=psi_W * psi_W'
rho_total = kron(rho_inF , rho_inF )

%% Display Non-Zero Diagonal Elements
diag_elems = diag(rho_total);
fprintf('Non-zero diagonal elements:\n');
for idx = 1:length(diag_elems)
    if ~isequal(diag_elems(idx), sym(0))
        fprintf('rho(%d,%d) = %s\n', idx, idx, char(diag_elems(idx)));
    end
end

Rc=CNOTn(6,1,4)*CNOTn(6,2,5)*CNOTn(6,3,6)*rho_total*(CNOTn(6,1,4)*CNOTn(6,2,5)*CNOTn(6,3,6))';

% Define computational basis projectors for sacrificial qubits
proj0 = [1, 0; 0, 0];  % |0⟩⟨0|
proj1 = [0, 0; 0, 1];  % |1⟩⟨1|


% Prepare basis combinations from 000 to 110
outcomes = dec2bin(0:6, 3);  % 000 to 110

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
    %disp(rho_norm);
    disp(rho_purified);
    
    fid=simplify(trace(rho_norm*rw));
    disp('Fidelity:');
    disp(fid);
   
end
