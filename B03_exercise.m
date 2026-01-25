% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------
% ----------------------------------------
% Build a model that "focuses" the wavefield in one direction and/or
% against one target.
% You can use array of sources  properly delayed (beam-forming), 
% and/or obstacles (reflectors) properly shaped (acoustic lens/horn)
%
% https://en.wikipedia.org/wiki/Architectural_acoustics

% -------------------------------------------------------------------------
%   Acoustic wave equation finite difference simulator
%   TASK: Build a model that focuses the wavefield on a target
%   Using: (1) source array + delays (beamforming)
%          (2) acoustic lens (velocity anomaly shaped)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  Acoustic wave equation finite difference simulator
%  TASK: Focus the wavefield on a target
%  Using: (1) source array + delays (beamforming)
%         (2) acoustic lens / horn (shaped velocity anomaly)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Acoustic wave equation finite difference simulator
% TASK: Focus the wavefield on a target using:
% (1) Source array with beamforming delays
% (2) Shaped velocity anomaly (acoustic lens / horn)
% -------------------------------------------------------------------------

clc; clear; close all;

% MODEL GRID

model.x = 0:2:2400;     % dx = 2 m
model.z = 0:2:1200;     % dz = 2 m
Nx = numel(model.x);
Nz = numel(model.z);
[X,Z] = meshgrid(model.x, model.z);

% BACKGROUND VELOCITY

v0 = 3000;
model.vel = v0 * ones(Nz, Nx);

% TARGET
target_x = 1600;
target_z = 700;

% ACOUSTIC LENS
lens_center_x = 1200;
lens_center_z = 500;
lens_a = 320;
lens_b = 220;

lens_mask = ((X - lens_center_x).^2)/(lens_a^2) + ...
            ((Z - lens_center_z).^2)/(lens_b^2) <= 1;

v_lens = 2000;
model.vel(lens_mask) = v_lens;

% SOURCE ARRAY + BEAMFORMING
source.x = 600:20:1100;
Ns = numel(source.x);
source.z = 80 * ones(1, Ns);

source.f0   = 18 * ones(1, Ns);
source.type = ones(1, Ns);
source.amp  = tukeywin(Ns, 0.4)';

v_ref = 2100;

t0_raw = zeros(1, Ns);
for i = 1:Ns
    dx = source.x(i) - target_x;
    dz = source.z(i) - target_z;
    dist = sqrt(dx^2 + dz^2);
    t0_raw(i) = dist / v_ref;
end
source.t0 = t0_raw - min(t0_raw);

% RECEIVERS
model.recx  = [target_x-300, target_x, target_x+300];
model.recz  = [target_z,     target_z, target_z];
model.dtrec = 0.001;

% SIMULATION PARAMETERS
simul.borderAlg  = 1;
simul.timeMax    = 1.0;
simul.printRatio = 8;
simul.higVal     = 0.4;
simul.lowVal     = 0.05;
simul.bkgVel     = 1;
simul.cmap       = 'turbo';

% RUN SIMULATION
recfield = acu2Dpro(model, source, simul);

% VELOCITY MODEL PLOT
figure;
imagesc(model.x, model.z, model.vel);
set(gca,'YDir','normal'); axis tight;
colormap(turbo); colorbar;
hold on;
plot(source.x, source.z, 'w.', 'MarkerSize', 10);
plot(model.recx, model.recz, 'ks', 'MarkerFaceColor','y');
plot(target_x, target_z, 'rp', 'MarkerSize', 14, 'MarkerFaceColor','r');
title('Velocity model, source array, receivers, and target');
xlabel('x (m)'); ylabel('z (m)');

% TRACE PLOT
t  = recfield.time;
tr = recfield.data;

figure;
plot(t, tr(:,1), 'LineWidth', 1.3); hold on;
plot(t, tr(:,2), 'LineWidth', 1.6);
plot(t, tr(:,3), 'LineWidth', 1.3);
grid on; axis tight;
xlabel('Time (s)');
ylabel('Amplitude');
title('Traces: before / at / after target');
legend('Before target','At target','After target','Location','best');

mid = 2;
[~, idx0] = max(abs(tr(:,mid)));

dt  = model.dtrec;
win = round(0.03 / dt);   % 60 ms window
i1 = max(1, idx0 - win);
i2 = min(length(t), idx0 + win);

Ewin_before = sum(tr(i1:i2,1).^2);
Ewin_at     = sum(tr(i1:i2,2).^2);
Ewin_after  = sum(tr(i1:i2,3).^2);

Gain_win   = Ewin_at / ((Ewin_before + Ewin_after)/2);
FocusRatio = Ewin_at / (Ewin_before + Ewin_at + Ewin_after);

% SANITY CHECKS
FullE = sum(tr.^2) * dt;
PeakA = max(abs(tr));

% PRINT RESULTS
fprintf('\n---- Windowed focus metric ----\n');
fprintf('v_ref = %.0f m/s\n', v_ref);
fprintf('Ewin before / at / after : %.3e  %.3e  %.3e\n', ...
        Ewin_before, Ewin_at, Ewin_after);
fprintf('Gain_win = %.2f | FocusRatio = %.2f\n', Gain_win, FocusRatio);

fprintf('\n---- Sanity checks ----\n');
fprintf('FullE before / at / after : %.3e  %.3e  %.3e\n', FullE);
fprintf('Peak before / at / after  : %.3e  %.3e  %.3e\n', PeakA);
fprintf('Peak ratio (at/before)    : %.2f\n', PeakA(2)/PeakA(1));
fprintf('FullE ratio (at/before)   : %.2f\n', FullE(2)/FullE(1));
