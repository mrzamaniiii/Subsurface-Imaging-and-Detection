% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------

% ----------------------------------------
% Practice with the examples (files named "A**_exercise.m"), 
% by changing source/receiver position, source type, velocity values

clear all;
clc;

% MODEL PARAMETERS

model.x = 0:1:1000;      % x axis (m)
model.z = 0:1:500;       % z axis (m)

Nx = numel(model.x);
Nz = numel(model.z);

% Velocity model
model.vel = zeros(Nz, Nx);

for kx = 1:Nx
    for kz = 1:Nz
        if model.z(kz) < 200
            model.vel(kz,kx) = 1500;   % layer 1
        else
            model.vel(kz,kx) = 2500;   % layer 2
        end
    end
end

% SOURCE PARAMETERS

source.x    = 400;     % source position
source.z    = 50;
source.f0   = 10;
source.t0   = 0.04;
source.amp  = 1;
source.type = 1;       % 1 = Ricker, 2 = sinusoidal

% --------------------------------------------------
% 3. RECEIVERS (change me!)
% --------------------------------------------------

model.recx  = 100:50:900;      % receiver line
model.recz  = ones(size(model.recx)) * 20;
model.dtrec = 0.004;

% --------------------------------------------------
% 4. SIMULATION PARAMETERS
% --------------------------------------------------

simul.borderAlg = 1;
simul.timeMax   = 0.5;
simul.printRatio= 10;
simul.higVal    = 0.3;
simul.lowVal    = 0.05;
simul.bkgVel    = 1;
simul.cmap      = 'gray';

% --------------------------------------------------
% 5. RUN SIMULATION
% --------------------------------------------------

recfield = acu2Dpro(model, source, simul);

% --------------------------------------------------
% 6. PLOT SEISMIC TRACES
% --------------------------------------------------

figure;
scal   = 1;
pltflg = 0;
scfact = 1;
colour = '';
clip   = [];

seisplot2(recfield.data, recfield.time, [], scal, pltflg, scfact, colour, clip);
xlabel('Receiver number');
ylabel('Time (s)');
title('Practice seismic traces');

% --------------------------------------------------
% PLOT VELOCITY MODEL + SOURCE + RECEIVERS
% --------------------------------------------------
figure;

imagesc(model.x, model.z, model.vel);
set(gca,'YDir','normal'); 
axis tight;

colormap(turbo);
colorbar;

hold on;

% Interface line (layer boundary at z = 200 m)
plot([model.x(1) model.x(end)], [200 200], 'k--', 'LineWidth', 1.2);

% Source
plot(source.x, source.z, 'rp', ...
     'MarkerFaceColor','r', 'MarkerSize', 10);

% Receivers
plot(model.recx, model.recz, 'w.', 'MarkerSize', 10);

xlabel('x (m)');
ylabel('z (m)');
title('Velocity model (two layers) with source and receivers');

legend({'Interface','Source','Receivers'}, ...
       'Location','northeast');
