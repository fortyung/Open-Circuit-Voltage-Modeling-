clear, clc
% Load data
load('C1209.mat');

% Initial data visualization
h = figure;
subplot(211); plot(Time, V, 'LineWidth', 2); grid on
xlabel('Time(hr)'); ylabel('Meas. Voltage(V)');
subplot(212); plot(Time, I, 'LineWidth', 2); grid on
xlabel('Time(hr)'); ylabel('Meas. Current(A)');
sgtitle('Batt: C1209');

% Question 1
% Calculate charge and discharge capacities

idxc = find(I >= 0);  % Indices for charge current
Ic = I(idxc);
Tc = Time(idxc);

idxd = find(I < 0);   % Indices for discharge current
Id = I(idxd);
Td = Time(idxd);

Qc = sum(diff(Tc) .* Ic(2:end)); % Charge capacity
Qd = -sum(diff(Td) .* Id(2:end)); % Discharge capacity

% ------------------------------------------

% Question 2: Compute SOC
SOC = zeros(length(Time), 1);
SOC(1) = 1; % Assume battery is full at start

for i = 2:length(I)
    if I(i) >= 0 
        SOC(i) = SOC(i-1) + (Time(i) - Time(i-1)) * (I(i)) / Qc;
    else
        SOC(i) = SOC(i-1) + (Time(i) - Time(i-1)) * (I(i)) / Qd;
    end
end

plot(Time, SOC)
grid on, axis tight
xlabel('Time (hr)')
ylabel('SOC')
title('Time VS SOC')


% ------------------------------------------

% Question 3: OCV model parameters

% Scaling factor
E = 0.175; 
zs = SOC * (1 - 2 * E) + E; % Scaled SOC

% First model: Combined model
P1 = [ones(length(V), 1), 1 ./ zs, zs, log(zs), log(1 - zs), I];
kest1 = (P1' * P1) \ (P1' * V); % LS estimate
R0est1 = kest1(end);
k1 = kest1(1:5);
OCV1 = k1(1) * ones(length(zs), 1) + k1(2) * (1 ./ zs) + k1(3) * zs + k1(4) * log(zs) + k1(5) * log(1 - zs);

% Second model: Combined+3 model
P = [ones(length(V), 1), 1 ./ zs, 1 ./ zs.^2, 1 ./ zs.^3, 1 ./ zs.^4, zs, log(zs), log(1 - zs), I];
kest = (P' * P) \ (P' * V); % LS estimate
R0est = kest(end);
k = kest(1:8);
OCV = k(1) * ones(length(zs), 1) + k(2) * (1 ./ zs) + k(3) * (1 ./ (zs.^2)) + ...
      k(4) * (1 ./ (zs.^3)) + k(5) * (1 ./ (zs.^4)) + k(6) * zs + k(7) * log(zs) + k(8) * log(1 - zs);

% Error Metrics Calculation
MEM1 = calculateErrorMetrics(V, OCV1, 6);
MEM2 = calculateErrorMetrics(V, OCV, 9);

% Display error metrics
disp('Error Metrics for Combined Model:');
disp(MEM1);
disp('Error Metrics for Combined+3 Model:');
disp(MEM2);

% Select best model based on error metrics
SOCx = 0:0.01:1;
zs = SOCx * (1 - 2 * E) + E; % Scaled SOC
zs = zs'; % Convert zs to column vector


% Calculate OCV for the SOC range
OCVx = k(1) * ones(length(zs), 1) + k(2) * (1 ./ zs) + k(3) * (1 ./ (zs .^ 2)) + ...
       k(4) * (1 ./ (zs .^ 3)) + k(5) * (1 ./ (zs .^ 4)) + k(6) * zs + ...
       k(7) * log(zs) + k(8) * log(1 - zs);

% Question 5

plot(zs,OCVx)
axis tight
xlabel('SOC')
ylabel('OCV in V')
title('OCV VS SOC')


% % Question 6
% (a 
function vio = socAtp(x,zs)
    ss = x * (1 - 2 * 0.175) + 0.175; % Scaled SOC
    vio = find(zs == ss); 
end

% OCV @ 25%
S25 = socAtp(0.25,zs);
OCV25 = OCVx(S25);

% OCV @50%
S50 = socAtp(0.50,zs);
OCV50 = OCVx(S50);

% OCV @ 75%
S75 = socAtp(0.75,zs);
OCV75 = OCVx(S75);

% (b)
% Find the SOC corresponding to the target voltage (4V)

OCV4 = min(OCVx(OCVx >= 4 & OCVx <5));
SOC4 = zs(OCVx == OCV4);

% Display the result
% % C)
v = 3.8;
e = v + 1*R0est;
target4 = e;
SOC3_8 = interp1(OCVx, zs, target4, "linear");


% Question 7: Open Circuit Voltage (OCV) Analysis Using the Combined Model

% Create properly sized SOC range for analysis
SOCx = (0:0.01:1)';  % Make sure it's a column vector
zs = SOCx * (1 - 2 * E) + E;

% Recalculate OCV1 for the new SOC range
P1_new = [ones(length(zs), 1), 1./zs, zs, log(zs), log(1-zs)];
OCV1_new = P1_new * k1;

% (a) Find OCV at specific SOC levels
ss25 = 0.25 * (1 - 2 * E) + E;
ss50 = 0.50 * (1 - 2 * E) + E;
ss75 = 0.75 * (1 - 2 * E) + E;

% Find closest indices
[~, idx25] = min(abs(zs - ss25));
[~, idx50] = min(abs(zs - ss50));
[~, idx75] = min(abs(zs - ss75));

% Get OCV values
OCV25_combined = OCV1_new(idx25);
OCV50_combined = OCV1_new(idx50);
OCV75_combined = OCV1_new(idx75);

% Display results
fprintf('Combined Model Results:\n');
fprintf('OCV @ 25%% SOC: %.4f V\n', OCV25_combined);
fprintf('OCV @ 50%% SOC: %.4f V\n', OCV50_combined);
fprintf('OCV @ 75%% SOC: %.4f V\n', OCV75_combined);

% (b) Find SOC corresponding to 4V
target_voltage = 4.0;
valid_indices = OCV1_new >= target_voltage;
if any(valid_indices)
    [~, min_idx] = min(abs(OCV1_new(valid_indices) - target_voltage));
    valid_positions = find(valid_indices);
    idx4V = valid_positions(min_idx);
    SOC4_combined = SOCx(idx4V);
    fprintf('SOC at 4V: %.4f\n', SOC4_combined);
else
    fprintf('No valid SOC found for 4V target\n');
end

% (c) Calculate SOC for charging to 3.8V
target_voltage_38 = 3.8;
target_ocv_combined = target_voltage_38 + R0est1;

% Use interpolation with the new arrays
SOC3_8_combined = interp1(OCV1_new, SOCx, target_ocv_combined, 'linear');
if ~isnan(SOC3_8_combined)
    fprintf('SOC for charging to 3.8V: %.4f\n', SOC3_8_combined);
else
    fprintf('Could not interpolate SOC for 3.8V target\n');
end




% Function for calculating error metrics
function metrics = calculateErrorMetrics(V, OCV, Nn)
    N = length(V);
    
    % R² fit
    V_mean = mean(V);
    R2 = (1 - (norm(OCV - V)^2 / norm(V - V_mean)^2)) * 100;
    
    % Max Error
    ME = max(abs(V - OCV));
    
    % Root Mean Square Error
    RMS = norm(V - OCV) / sqrt(N - Nn);
    
    % Store results in a struct
    metrics = struct('R2', round(R2, 4), 'ME', round(ME, 4), 'RMS', round(RMS, 4));
end