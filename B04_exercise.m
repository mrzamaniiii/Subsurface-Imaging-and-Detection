% ----------------------------------------
% Refraction  (with large borders to keep noise away)
% ----------------------------------------
%
% -------------------------------------------------------------------------
%       Acoustic wave equation finite diference simulator
% -------------------------------------------------------------------------
%
% This file is a manual and example script using the acoustic simulator
% The user defines the velocity model, the source parameters and 
% the simulation parameters
% The program displays the "evolving" acoustic pressure field amplitude
% in a movie-like figure (snapshots)
% The user can define an optional set of receivers positions: in this case
% the program shows and outputs the pressure recorded at the receivers
% (seismic traces). The receivers that are out of the model are
% automatically rejected, and so the seismic traces can be less than the
% input receivers: the output structure contains the actual position of the
% 'valid' receivers
%
% -------------------------------------------------------------------------
% 1. Model parameters in structure 'model'
%
% Compulsory parameters
% model.x: vector of x grid coordinates [m], Nx elements
% model.z: vector of z grid coordinates [m], Nz elements
% model.vel: matrix of velocity values [m/s], (Nz,Nx) elements
%
% Optional parameters
% model.recx: vector of x coordinates of receivers [m], Nr elements
% model.recz: vector of z coordinates of receivers [m], Nr elements
% model.dtrec: max time sampling interval for seismic traces [s]
%
% -------------------------------------------------------------------------
% 2. Source parameters in structure 'source'
% source coordinates are rounded to nearest grid point
% all of them can be **vectors**, in order to simulate multiple sources
%
% Compulsory parameters
% source.x:  x coordinate of source [m]
% source.z:  z coordinate of source [m]
% source.f0: central frequency of source Ricker wavelet [Hz]
% source.t0: time of source emission (referred to max peak of Ricker)
% source.type: 1 is Ricker, 2 is sinusoid at frequency source.f0
% source.amp: multiplier of source amplitude
%
% -------------------------------------------------------------------------
% 3. Simulation and graphic parameters in structure 'simul'
%
% simul.timeMax:     max simulation time [s]
% simul.borderAlg:   Absorbing boundaries (Yes:1, No:0)
% simul.printRatio:  pressure map shown every printRatio comput. time steps
% simul.higVal:  colormap between -highVal and + highVal (from 0 to 1)
% simul.lowVal:  values between -lowVal and +lowVal zeroed (from 0 to 1)
% simul.bkgVel:  velocity matrix as  a "shadow" in the images (1:yes, 0:no)
%
% Optional parameters
% simul.cmap:    colormap (default 'gray')
%
% -------------------------------------------------------------------------
% 4. Acoustic simulator program call
%
%  recfield=acu2Dpro(model,source,simul);
%
%  recfield.time: time axis of recorded signal [s], Nt elements
%  recfield.data: matrix of pressure at the receivers, (Nt,Nr1)
%  recfield.recx: vector of x grid coordinates of receivers [m], Nr1 elements
%  recfield.recz: vector of z grid coordinates of receivers [m], Nr1 elements
%
% -------------------------------------------------------------------------
% 5. Plotting seismic trace program call
%
%  seisplot2(recfield.data,recfield.time)
% 
%   [fact]=seisplot2(datain,t,tr,scal,pltflg,scfact,colour,clip)
%  
%   function for plotting seismic traces
%  
%   INPUT
%   datain  - input matrix of seismic traces
%   t       - time axis
%   tr      - trace axis
%   scal    - 1 for global max, 0 for global ave, 2 for trace max
%   pltflg  - 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
%             2 plots wiggle traces only, 3 imagesc gray, 4 pcolor gray
%   scfact  - scaling factor
%   colour  - trace colour, default is black
%   clip    - clipping of amplitudes (if <1); default no clipping
%  
%   OUPTPUT
%   fact    - factor that matrix was scaled by for plotting
%   if you want to plot several matrices using the same scaling factor,
%   capture 'fact' and pass it as input variable 'scal' to future matrices
%   with 'scfact' set to 1
% 
% -------------------------------------------------------------------------
% 6. Essential theory
%
%  Sampling interval for stability (computed automatically)
%  dt=0.95*sqrt(1/2)*min(dx,dz)/vMax
%
%  Courant-Friedrick-Levy (CFL) stability condition (display warning)
%  dx<(vMin/(f0*2*8)))
%  dz<(vMin/(f0*2*8)))
%
%  -----------
%  Ricker wavelet
%
%                  *                        !    ***
%                 * *                       !   * ! *
%    ____________*___*___________ T         !  *  !   *
%    ********   *     *   *******           ! *   !     *
%            ***       ***                  !*    !        * * * *
%             >   TD    <                   ------+---------------- F
%                                                F0
%
%    s(t) = (1-Y2*T*T)*exp(-Y2*T*T/2)   S(f) = 2*F^2/(f0^3*sqrt(pi))*
%    Y2 = 2*pi^2*f0^2                     *exp(-F*F/(f0^2))
%    TD = sqrt(6)/(pi*f0)
%
%  -----------


% -------------------------------------------------------------------------
% Refraction example (with large borders to keep boundary noise away)
% Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Refraction example (with large borders to keep boundary noise away)
% Acoustic wave equation finite difference simulator
% Fixes applied:
%   (1) Shot gather x-axis = offset (m), not receiver number
%   (2) Symmetric receiver spread around the source (two-sided offsets)
% -------------------------------------------------------------------------

clc; clear; close all;

dx = 2; dz = 2;

x_main_min = 0;   x_main_max = 2000;
z_main_min = 0;   z_main_max = 1000;

x_border = 600;
z_border = 400;

model.x = (x_main_min - x_border):dx:(x_main_max + x_border);
model.z = (z_main_min):dz:(z_main_max + z_border);

[X,Z] = meshgrid(model.x, model.z);

v1 = 1500;
v2 = 3000;

z0  = 350;
dip = 0.08;
x0  = (x_main_min + x_main_max)/2;

z_interface = z0 + dip*(X - x0);

model.vel = v1 * ones(size(Z));
model.vel(Z >= z_interface) = v2;

source.x    = 1000;
source.z    = 60;
source.f0   = 12;
source.t0   = 0.05;
source.type = 1;
source.amp  = 1;

rec_dx  = 20;
max_off = 900;
offsets = (-max_off:rec_dx:max_off);

model.recx  = source.x + offsets;
model.recz  = source.z * ones(size(model.recx));
model.dtrec = 0.002;

simul.timeMax    = 1.2;
simul.borderAlg  = 1;
simul.printRatio = 10;
simul.higVal     = 0.35;
simul.lowVal     = 0.03;
simul.bkgVel     = 1;
simul.cmap       = 'gray';

recfield = acu2Dpro(model, source, simul);

figure;
imagesc(model.x, model.z, model.vel);
set(gca,'YDir','normal'); axis tight;
colormap(turbo); colorbar; hold on;

hS = plot(source.x, source.z, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 10);
hR = plot(model.recx, model.recz, 'w.', 'MarkerSize', 8);

x_line = model.x;
z_line = z0 + dip*(x_line - x0);
plot(x_line, z_line, 'k--', 'LineWidth', 1.2);

xlabel('x (m)'); ylabel('z (m)');
title('Refraction model + source/receivers');
legend([hS hR], {'Source','Receivers'}, 'Location','northeast');

offset_valid = recfield.recx - source.x;

% Sort by offset so the x-axis is always monotonic
[off_sorted, idx] = sort(offset_valid);
data_sorted = recfield.data(:, idx);

figure;
scal   = 1;
pltflg = 0;
scfact = 1;
colour = '';
clip   = 0.9;

seisplot2(data_sorted, recfield.time, off_sorted, scal, pltflg, scfact, colour, clip);
xlabel('Offset (m)');
ylabel('Time (s)');
title('Shot gather (refraction): direct + refracted arrivals');
grid on;
