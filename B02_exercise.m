% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------

% --------------------------------------------------
% In a homogeneous medium, compute the energy of traces at increasing 
% distance from the source
% Plot the resulting energy vs distance

clear all;
clc;

% MODEL PARAMETERS (Homogeneous medium)

model.x = 0:1:1200;     
model.z = 0:1:600;      

Nx = numel(model.x);
Nz = numel(model.z);

model.vel = 2000 * ones(Nz, Nx);   % homogeneous velocity

% SOURCE PARAMETERS

source.x    = 600;     % source in the center
source.z    = 100;
source.f0   = 20;      % Hz
source.t0   = 0.04;
source.amp  = 1;
source.type = 1;       % Ricker

% --------------------------------------------------
% 3. RECEIVERS (increasing distance from source)
% --------------------------------------------------

model.recx  = 650:50:1150;   % increasing offsets
model.recz  = ones(size(model.recx)) * 100;
model.dtrec = 0.002;

offsets = abs(model.recx - source.x);  % distance from source

% --------------------------------------------------
% 4. SIMULATION PARAMETERS
% --------------------------------------------------

simul.borderAlg = 1;
simul.timeMax   = 0.6;
simul.printRatio = 10;
simul.higVal = 0.3;
simul.lowVal = 0.05;
simul.bkgVel = 1;
simul.cmap = 'gray';

% --------------------------------------------------
% 5. RUN SIMULATION
% --------------------------------------------------

recfield = acu2Dpro(model, source, simul);

% --------------------------------------------------
% 6. ENERGY COMPUTATION
% --------------------------------------------------

dt = recfield.time(2) - recfield.time(1);
Nr = size(recfield.data, 2);

energy = zeros(1, Nr);

for i = 1:Nr
    trace = recfield.data(:, i);
    energy(i) = sum(trace.^2) * dt;   % energy of trace
end

% --------------------------------------------------
% 7. PLOTS
% --------------------------------------------------

% Seismic traces
figure;
seisplot2(recfield.data, recfield.time, offsets, 2, 0, 1);
xlabel('Offset (m)');
ylabel('Time (s)');
title('Seismic Traces at Increasing Distance');

% Energy vs Distance
figure;
plot(offsets, energy, 'o-', 'LineWidth', 2);
xlabel('Distance from Source (m)');
ylabel('Energy');
title('Energy Decay with Distance (Homogeneous Medium)');
grid on;
