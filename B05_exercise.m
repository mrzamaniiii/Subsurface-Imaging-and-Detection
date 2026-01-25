% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------
% Create your own model and example....
% --------------------------------------------------

% -------------------------------------------------------------------------
% Acoustic wave equation finite difference simulator
% Example: refraction + reflection/scattering from an obstacle
% Fix applied:
%   (1) Put source/receivers near the surface (avoid bottom boundary effects)
% -------------------------------------------------------------------------

clc; clear; close all;

% Grid
model.x = 0:1:1000;
model.z = 0:1:500;

Nx = numel(model.x);
Nz = numel(model.z);

[X,Z] = meshgrid(model.x, model.z);

% Velocity model
model.vel = zeros(Nz, Nx);

% Layer 1: z < 150
model.vel(Z < 150) = 800;

% Layer 2: 150 <= z < 350 : v = 1000 + 2*(z-150)
mask2 = (Z >= 150) & (Z < 350);
model.vel(mask2) = 1000 + 2*(Z(mask2) - 150);

% Layer 3: z >= 350
model.vel(Z >= 350) = 2200;

% High-speed obstacle (reflection + scattering)
model.vel(250:300, 400:600) = 4000;

% Source (moved near surface)
source.x = 500;
source.z = 60;
source.f0 = 10;
source.t0 = 0.05;
source.amp = 1;
source.type = 1;

% Receivers (moved near surface)
model.recx = 100:50:900;
model.recz = 60 * ones(size(model.recx));
model.dtrec = 0.002;

% Simulation
simul.borderAlg = 1;
simul.timeMax = 0.6;
simul.printRatio = 10;
simul.higVal = 0.3;
simul.lowVal = 0.02;
simul.bkgVel = 1;
simul.cmap = 'gray';

% Run
recfield = acu2Dpro(model, source, simul);

% Plot velocity model + geometry
figure;
imagesc(model.x, model.z, model.vel);
set(gca,'YDir','normal'); axis tight;
colormap(turbo); colorbar;
hold on;
hS = plot(source.x, source.z, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 10);
hR = plot(model.recx, model.recz, 'w.', 'MarkerSize', 10);
xlabel('x (m)'); ylabel('z (m)');
title('Velocity model + source + receivers');
legend([hS hR], {'Source','Receivers'}, 'Location','northeast');

% Shot gather (unchanged: receiver number on x-axis)
figure;
scal = 1;
pltflg = 0;
scfact = 1;
colour = '';
clip = 0.9;

seisplot2(recfield.data, recfield.time, [], scal, pltflg, scfact, colour, clip);
xlabel('Receiver Number');
ylabel('Time (s)');
title('Shot gather: refraction + reflection/scattering from obstacle');
