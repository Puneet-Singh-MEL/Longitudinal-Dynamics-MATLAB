% Main script for Aircraft Longitudinal Dynamics Analysis

% Clear workspace, close figures, and clear command window
clear all; 
close all; 
clc;

%% Section 2.1 Experimental Data
% Import Data
dataFiles = {'Step_Output.txt', 'Step_Input.txt', 'Singlet_Output.txt', ...
             'Singlet_Input.txt', 'Doublet_Output.txt', 'Doublet_Input.txt', ...
             'Slowramp_Output.txt', 'Slowramp_Input.txt', 'Fastramp_Output.txt', ...
             'Fastramp_Input.txt', 'Arb_Output.txt', 'Arb_Input.txt'};
data = importDataFiles(dataFiles);

% Plot Data
plotData(data);

%% Section 3.1.1 Stability Derivatives
% Calculate stability derivatives based on technical characteristics
stabilityParams = calculateStabilityDerivatives();

%% Section 3.1.2 Control Derivative
% Calculate control derivatives
controlParams = calculateControlDerivatives();

%% Section 3.2.1 Eigenvalues
% Calculate eigenvalues of the system
eigenvalues = calculateEigenvalues(stabilityParams);

%% Section 3.3 Aircraft Response Transfer Function
% Define Transfer Function
transferFunction = defineTransferFunction();

%% Section 3.5.1 Step Input
% Plot Theoretical and Experimental Step Input
plotStepInput();

%% Function Definitions

% Function to import data from specified files
function data = importDataFiles(dataFiles)
    % Initialize a structure to hold data
    data = struct();
    
    % Loop through each file and import data
    for i = 1:length(dataFiles)
        [~, name, ~] = fileparts(dataFiles{i});
        try
            data.(name) = importdata(dataFiles{i});
        catch ME
            fprintf('Error importing file %s: %s\n', dataFiles{i}, ME.message);
            % Handle the error gracefully or rethrow it as needed
            rethrow(ME);
        end
    end
end

% Function to plot input and output data
function plotData(data)
    % Loop through each signal and plot input and output
    signals = fieldnames(data);
    for i = 1:2:length(signals)
        figure;
        % Plot raw data
        plot(linspace(0, max(data.(signals{i}).data(:,1)), length(data.(signals{i+1}).data)), data.(signals{i+1}).data);
        hold on;
        % Filter and plot filtered data
        filteredSig = filterSignal(data.(signals{i}).data, 5, 10000);
        title(strrep(signals{i}, '_', ' '));
    end
end

% Function to filter a signal using a bandpass filter
function filteredSig = filterSignal(raw_signal, fmax, samprate)
    try
        % Extract time and pitch data
        t = raw_signal(:,1);
        pitch = raw_signal(:,2);
        
        % Design and apply bandpass filter
        cut_off = [0.01 fmax];
        d = designfilt('bandpassiir', 'FilterOrder', 2, 'HalfPowerFrequency1', cut_off(1), 'HalfPowerFrequency2', cut_off(2), 'SampleRate', samprate);
        filteredSig = filter(d, pitch);
        
        % Plot original and filtered signals
        plot(t, pitch);
        hold on;
        plot(t, filteredSig, 'r', 'LineWidth', 2);
        legend('Input', 'Raw Output', 'Filtered Output');
        set(gca, 'FontSize', 18, 'FontName', 'Times New Roman');
        xlabel('Time', 'FontName', 'Times New Roman');
        ylabel('Pitch Response', 'FontName', 'Times New Roman');
    catch ME
        fprintf('Error filtering signal: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end

% Function to calculate stability derivatives
function stabilityParams = calculateStabilityDerivatives()
    try
        % Declare technical characteristics and parameters
        C_l_a_wing = 0.108/(pi/180); 
        C_l_a_tail = 0.108/(pi/180); 
        AR_wing = 4.544;
        AR_tail = 3.4; 
        c_r = 125/1000; 
        c_t = 94/1000; 
        l_tail = 339/1000; 
        S_tail = 0.014; 
        S_wing = 0.05425; 
        X_CoM = 10/1000; 
        X_ac = 27.3/1000; 
        n = 0.9; 
        de_da = 0.606; 
        u = 27; 
        Q = 0.5 * 1.225 * u^2; 
        I_yy = 0.039; 
        
        % Calculate derived parameters
        C_L_a_wing = C_l_a_wing / (1 + C_l_a_wing / (pi * AR_wing)); 
        C_L_a_tail = C_l_a_tail / (1 + C_l_a_tail / (pi * AR_tail)); 
        lamda = c_t / c_r; 
        c = c_r * (2 / 3) * ((1 + lamda + lamda^2) / (1 + lamda)); 
        V_H = (l_tail * S_tail) / (S_wing * c); 
        C_m_a = C_L_a_wing * ((X_CoM / c) - (X_ac / c)) - n * V_H * C_L_a_tail * (1 - de_da); 
        M_w = C_m_a * Q * S_wing * c / (u * I_yy); 
        M_a = u * M_w; 
        C_M_q = -2 * n * C_L_a_tail * V_H * l_tail / c; 
        M_q = C_M_q * (c / (2 * u)) * Q * S_wing * c / I_yy;
        
        % Store stability parameters in a structure
        stabilityParams = struct('M_w', M_w, 'M_a', M_a, 'M_q', M_q);
    catch ME
        fprintf('Error calculating stability derivatives: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end

% Function to calculate control derivatives
function controlParams = calculateControlDerivatives()
    try
        % Declare parameters
        u = 27; 
        Q = 0.5 * 1.225 * u^2; 
        t = 0.68; 
        C_l_a_tail = 0.108 / (pi / 180); 
        AR_tail = 3.4; 
        n = 0.9; 
        l_tail = 339 / 1000; 
        S_tail = 0.014; 
        S_wing = 0.05425; 
        c_r = 125 / 1000; 
        c_t = 94 / 1000; 
        lamda = c_t / c_r; 
        c = c_r * (2 / 3) * ((1 + lamda + lamda^2) / (1 + lamda)); 
        I_yy = 0.039; 
        C_L_a_tail = C_l_a_tail / (1 + C_l_a_tail / (pi * AR_tail)); 
        dCLTail_dDe = t * C_L_a_tail;
        V_H = (l_tail * S_tail) / (S_wing * c); 
        C_M_d_E = -n * V_H * dCLTail_dDe;
        M_d_E = C_M_d_E * (Q * S_wing * c) / I_yy;
        
        % Store control parameters in a structure
        controlParams = struct('M_d_E', M_d_E);
    catch ME
        fprintf('Error calculating control derivatives: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end

% Function to calculate eigenvalues of the system
function eigenvalues = calculateEigenvalues(stabilityParams)
    try
        % Calculate eigenvalues using symbolic math toolbox
        syms s 
        M_a = stabilityParams.M_a; 
        M_q = stabilityParams.M_q; 
        A = [0 1; M_a M_q]; 
        I = eye(2, 2); 
        T = det(I * s - A); 
        s = solve(T == 0); 
        
        % Create transfer function and plot pole-zero map
        G = tf([0], [1 -M_q -M_a]); 
        pzplot(G)
        
        % Store eigenvalues
        eigenvalues = s;
    catch ME
        fprintf('Error calculating eigenvalues: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end

% Function to define transfer function for aircraft response
function transferFunction = defineTransferFunction()
    try
        % Define transfer function parameters
        syms s;
        mq = -2.4058;
        ma = -121.8784;
        mde = -130.2904;

        % Define transfer function matrix A
        A = [s - mq 1; ma s] / (s^2 - mq * s - ma);

        % Define input matrix B and output matrix C
        B = [0; mde];
        C = [1 0];
        D = 0;

        % Compute transfer function G(s)
        G = C * inv(s * eye(size(A)) - A) * B + D;

        % Return symbolic transfer function
        transferFunction = G;
    catch ME
        fprintf('Error defining transfer function: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end

function plotStepInput()
    try
        % Load experimental step input data
        load('Step_Output.txt');
        Input = Step_Output(:,2);
        Time = Step_Output(:,1);

        % Plot experimental step input
        figure;
        plot(Time, Input);
        hold on;

        % Define and plot theoretical step input
        syms t;
        F(t) = (0.0913672841821045 * (t - 0) + (-0.055031)) * (heaviside(t - 0) - heaviside(t - 0.55)) + 0.0913672841821045 * (-0.055031) * (heaviside(t - 0.55) - heaviside(t - 1));
        fplot(F, [0 5], 'r', 'LineWidth', 2);

        % Add legend and labels
        legend('Experimental Input', 'Theoretical Input');
        xlabel('Time');
        ylabel('Input');
        hold off;
    catch ME
        fprintf('Error plotting step input: %s\n', ME.message);
        % Handle the error gracefully or rethrow it as needed
        rethrow(ME);
    end
end