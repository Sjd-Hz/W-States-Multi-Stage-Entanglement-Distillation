 %% (psi_1,psi_2) 
 % Transforming the GHZ-like state: |ψ⟩ = (|001⟩ + |010⟩)/√2 to
 %  W-like state: |ψ⟩ = (|001⟩ + |010⟩+ |100⟩+ |111⟩)/2
 %  via CNOT teleportation with 
 %\psi_2: bell between qubits 1 and 3, so CNOT between 1 and 3
    
    %% Initialization
    clc;
    clear;
    close all;
    
    %% Fundamental Definitions
    % Qubit state vectors
    ket0 = [1; 0];  % |0⟩
    ket1 = [0; 1];  % |1⟩
    
    % Identity operators
    Id2 = eye(2);
    Id4 = eye(4);
    Id8 = eye(8);
    
    % Single-qubit gates
    Had = (1/sqrt(2)) * [1 1; 1 -1];  % Hadamard gate
    X = [0 1; 1 0];                   % Pauli X gate
    Z = [1 0; 0 -1];                  % Pauli Z gate
    XZ = X * Z;                       % X followed by Z
    
    % Projection operators
    proj_p = 0.5 * [1 1; 1 1];        % |+⟩⟨+|
    proj_m = 0.5 * [1 -1; -1 1];      % |-⟩⟨-|
    proj_0 = [1 0; 0 0];              % |0⟩⟨0|
    proj_1 = [0 0; 0 1];              % |1⟩⟨1|

    %% State Preparation
    % Initial state (|001⟩ + |010⟩)/√2
    psi1 = (1/sqrt(2)) * (kron(ket0, kron(ket0, ket1)) + ...
            (kron(ket0, kron(ket1, ket0))));
    rho1 = psi1 * psi1'
    
    % Apply Hadamard to qubit 1
    Had1 = kron(Had,kron(Id2, Id2));
    rho_had = Had1 * rho1 * Had1';
    
    % Create Bell state (|01⟩ + |10⟩)/√2 between qubits B1 and B3
    bell = (1/sqrt(2)) * (kron(ket0, ket1) + kron(ket1, ket0));
    rho_bell = bell * bell';
    
    % Combine systems (now 5 qubits: 1,2,3,B1,B3)
    rho_total = kron(rho_had, rho_bell);

    %% CNOT Operations
    % CNOT: qubit 1 controls B1 (qubit 4)
    CNOT_1to4 = CNOTn(5,1,4);
    
    % CNOT: B3 (qubit 5) controls qubit 3
    CNOT_5to3 = CNOTn(5,5,3);
    
    % Apply both CNOTs
    rho_after_CNOTs = CNOT_5to3 * CNOT_1to4 * rho_total * CNOT_1to4' * CNOT_5to3';

    %% Measurement and Correction
    % Define measurement projectors for qubits 4 and 5 (B1 and B3)
    projectors = {
        '1+', kron(Id8, kron(proj_1, proj_p)), ...
        '1-', kron(Id8, kron(proj_1, proj_m)), ...
        '0+', kron(Id8, kron(proj_0, proj_p)), ...
        '0-', kron(Id8, kron(proj_0, proj_m))
    };

    % Correction operators for qubits 1, 2, 3
    corrections = {
        '1+', kron(kron(Id2, Id2), Id2), ...
        '1-', kron(kron(Z, Id2), Id2), ...
        '0+', kron(kron(Id2, Id2), X), ...
        '0-', kron(kron(Z, Id2), X)
   };

%% Step 4: Loop Over All Measurement Outcomes and Apply Corrections
for k = 1:2:length(corrections)
    label = corrections{k};
    corr_op = corrections{k+1};
    proj_op = projectors{k+1};

    % Apply measurement projector
    rho_proj = proj_op * rho_after_CNOTs * proj_op';

    % Trace out measured qubits (qubits 4 and 5)
    rho_reduced = PartialTrace(rho_proj, [4, 5]);

    % Apply correction operator
    rho_corrected = corr_op * rho_reduced * corr_op';

    % Normalize the final state
    prob_success = trace(rho_corrected);
    rho_final = rho_corrected / prob_success;

    % Display result
    fprintf('\n--- Measurement Outcome: |%s⟩ ---\n', label);
    fprintf('Success Probability: %.4f\n', double(prob_success));
    disp('Corrected and Normalized State ρ_final:');
    disp((rho_final));
end
