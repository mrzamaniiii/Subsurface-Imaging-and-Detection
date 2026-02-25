function recfield=acu2Dpro(model,source,simul);

%  Version 6 - Sep 2020
%  Version 7 - Mar 2021 (random signal) source
% -------------------------------------------------------------------------
%       Acoustic wave equation finite diference simulator
% -------------------------------------------------------------------------
%
% The user defines the velocity model, the source and simulation parameters
% The program displays the "evolving" acoustic pressure field amplitude
% in a movie-like figure (snapshots)
% The user can define an optional set of receivers positions: in this case
% the program outputs the pressure recorded at the receivers
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
% model.recx: vector of x grid coordinates of receivers [m], Nr elements
% model.recz: vector of z grid coordinates of receivers [m], Nr elements
% model.dtrec: max time sampling for rec traces [s]
%
% -------------------------------------------------------------------------
% 2. Source parameters in structure 'source'
% source coordinates are rounded to nearest grid point
% all of them can be vectors, in order to simulate multiple sources
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
%-------------%

[Nz,Nx] = size(model.vel);      %calcolo le dimensioni del modello

m0   = zeros(Nz,Nx);      % matrice di simulazione
m1   = zeros(Nz,Nx);      % matrice di simulazione
m2   = zeros(Nz,Nx);      % matrice di simulazione

model.Nx = Nx;
model.Nz = Nz;

vMax = max(max(model.vel));           % velocità massima
vMin = min(min(model.vel));           % velocità minima

dx = model.x(2)-model.x(1);
x0 = model.x(1);
dz = model.z(2)-model.z(1);
z0 = model.z(1);

f0 = source.f0;
t0 = source.t0;
Nsource = numel(source.amp);


dt=0.95*sqrt(1/2)*min(dx,dz)/vMax;  % Campionamento temporale, dalla formula di stabilità

% round source position to nearest point in the velocity grid
for kkk=1:Nsource,
  [MM,Ix_source(kkk)]=min(abs(source.x(kkk)-model.x));
  [MM,Iz_source(kkk)]=min(abs(source.z(kkk)-model.z));
end

recfield               = {};
model.compute_recfield = 0;
saveratio              = 0;
Nr                     = 0;
Nr1                    = 0;

if isfield(model,'recx'),
  Nr  = numel(model.recx);
  Nr1 = 0;
  for k=1:Nr;
    
    if (model.recx(k)>max(model.x)) || (model.recx(k)<min(model.x)) || ...
        (model.recz(k)>max(model.z)) || (model.recz(k)<min(model.z)),
      disp('Receiver out of model removed')
      continue
    end
    
    %     model.recIx(k)=Ix_rec;
    %     model.recIz(k)=Iz_rec;
    Nr1=Nr1+1;
    [MM,Ix_rec(Nr1)]=min(abs(model.recx(k)-model.x));
    [MM,Iz_rec(Nr1)]=min(abs(model.recz(k)-model.z));
    
    recfield.recx(Nr1)=model.x(Ix_rec(Nr1));
    recfield.recz(Nr1)=model.z(Iz_rec(Nr1));
    
  end
  
  if Nr1>0;
    model.compute_recfield=1;
    
    saveratio=max([1 floor(model.dtrec/dt)]);
    Ntrec=round(simul.timeMax/dt/saveratio);
    
    recfield.data=zeros(Ntrec,Nr1);
    recfield.time=zeros(1,Ntrec);
  end
end

if isfield(simul,'cmap')==0,
  simul.cmap='gray';
end

%----------------
%-- Check data --
%----------------
if  ((dx>(vMin/(max(f0)*2*8)))||(dz>(vMin/(max(f0)*2*8))))
  disp('************** WARNING: CFL condition not fulfilled');
end

for kkk=1:Nsource,
  if (source.x(kkk)>=max(model.x)) || (source.x(kkk)<=min(model.x)) || ...
      (source.z(kkk)>=max(model.z)) || (source.z(kkk)<=min(model.z))
    disp('***** Error: source on the border or out of model');
    return
  end
end

Idx_t  = 0;                   %Idx_t-esimo istante della simulazione
t_rec  = 0;
t_curr = 0;                  %tempo parziale della simulazione [s]

for kkk=1:Nsource,
  if source.type(kkk)==2 | source.type(kkk)==3 ,
    t_source_end(kkk)=simul.timeMax;
    t_source_start(kkk)=t0(kkk);
  elseif source.type(kkk)==1,
    t_source_start(kkk)=t0(kkk)-2*sqrt(6)/(pi*f0(kkk));
    t_source_end(kkk) = 2*sqrt(6)/(pi*f0(kkk))+t0(kkk);
  end
end

%--------------------------------------------------------------------

nv=zeros(size(model.vel));
if (simul.bkgVel==1)
  nv=-(model.vel/max(max(model.vel)))*simul.higVal/3;
end;

figure;
hh.ax=axes('units','normalized','Position',[0.1 0.2 0.8 0.7],'fontsize',8);
hh.text=uicontrol('style','text','units','normalized',...
  'Position',[0.1 0.01 .3 .06],...
  'string','');
hh.pause=uicontrol('style','togglebutton','units','normalized',...
  'Position',[0.45 0.01 .2 .08],...
  'value',0,'string','Pause');
hh.end=uicontrol('style','togglebutton','units','normalized',...
  'Position',[0.7 0.01 .2 .08],...
  'value',0,'string','End');


hh.img=imagesc(model.x,model.z,zeros(size(model.vel)),'parent',hh.ax);
colormap(hh.ax,simul.cmap);

caxis(hh.ax,[-simul.higVal simul.higVal]);
axis(hh.ax,'image');
xlabel(hh.ax,'horizontal axis [m]')
ylabel(hh.ax,'vertical axis [m]')
hold(hh.ax,'on')
if model.compute_recfield==1,
  plot(hh.ax,recfield.recx,recfield.recz,'.r','markersize',8)
end
set(hh.ax,'fontsize',8)
plot(hh.ax,source.x,source.z,'*r','markersize',8)

%--------
%- Core -
%--------
while (t_curr<=simul.timeMax)
  
  
  while hh.pause.Value==1;
    pause(.4)
  end
  
  if hh.end.Value==1;
    break
  end
  
  
  %Calcolo l'ondina di Ricker all'istante corrente
  
  for kkk=1:Nsource,
    app1=Idx_t.*dt-t0(kkk);
    
    if source.type(kkk)==1,
      y1=2.*(pi.*f0(kkk)).^2;
      r1(kkk)=(1.-y1.*(app1.^2)) .* exp(-y1.*(app1.^2)./2) ;
    elseif source.type(kkk)==2,
      r1(kkk)=sin(2*pi*f0(kkk)*app1);
    elseif source.type(kkk)==3,
      r1(kkk)=randn(1,1);
    end
    
    r1(kkk)=r1(kkk)*source.amp(kkk);
  end
  
  %Calcolo i valori della matrice e applico l'ondina di Ricker
  if(mod(Idx_t,3)==0)
    [m0 m1 m2]=fde(m0,m1,m2,model.vel,Nx,Nz,dz,dx,dt,simul.borderAlg);
    for kkk=1:Nsource,
      if t_curr>t_source_start(kkk) & t_curr<t_source_end(kkk),
        
        m2(Iz_source(kkk),Ix_source(kkk))=r1(kkk); %m2(Iz_source,Ix_source)+r1;
      end
    end
    
  elseif(mod(Idx_t,3)==1)
    [m1 m2 m0]=fde(m1,m2,m0,model.vel,Nx,Nz,dz,dx,dt,simul.borderAlg);
    for kkk=1:Nsource,
      if t_curr>t_source_start(kkk) & t_curr<t_source_end(kkk),
        
        m0(Iz_source(kkk),Ix_source(kkk))=r1(kkk); %m2(Iz_source,Ix_source)+r1;
      end
    end
  elseif(mod(Idx_t,3)==2)
    [m2 m0 m1]=fde(m2,m0,m1,model.vel,Nx,Nz,dz,dx,dt,simul.borderAlg);
    for kkk=1:Nsource,
      if t_curr>t_source_start(kkk) & t_curr<t_source_end(kkk),
        m1(Iz_source(kkk),Ix_source(kkk))=r1(kkk); %m2(Iz_source,Ix_source)+r1;
      end
    end
  end
  % figure(1000);
  % plot(t_curr,m1(Ix_source,Iz_source),'o',t_curr,r1,'*')
  % figure(1)
  %Stampo
  
  % -------------------------------------------------------------------
  % plot snapshot
  
  if (mod(Idx_t,simul.printRatio)==0)
    if(mod(Idx_t,3)==0)
      set(hh.img,'CData',m2+nv);
    elseif(mod(Idx_t,3)==1)
      set(hh.img,'CData',m0+nv);
    elseif(mod(Idx_t,3)==2)
      set(hh.img,'CData',m1+nv);
    end
    
    
    timeParStr=sprintf( 'Time = %0.5f [s] of %0.1f', [t_curr simul.timeMax]);
    hh.text.String=timeParStr;
    pause(1e-4)
  end
  drawnow
  % -------------------------------------------------------------------
  % save receiver
  
  %   if model.compute_recfield==1 & (mod(Idx_t,saveratio)==0) ,
  %     if(mod(Idx_t,3)==0)
  %       recsample=m2(sub2ind([Nz Nx],Iz_rec,Ix_rec));
  %     elseif(mod(Idx_t,3)==1)
  %       recsample=m0(sub2ind([Nz Nx],Iz_rec,Ix_rec));
  %     elseif(mod(Idx_t,3)==2)
  %       recsample=m1(sub2ind([Nz Nx],Iz_rec,Ix_rec));
  %     end
  %
  %     t_rec=t_rec+1;
  %     recfield.data(t_rec,:) = recsample(:)';
  %     recfield.time(t_rec)   = t_curr(:);
  %
  %     timeParStr=sprintf( 'Time = %0.5f [s]', t_curr);
  %   hh.text.String=timeParStr;
  %   pause(1e-4)
  %   drawnow
  %   end
  
  % -------------------------------------------------------------------
  % new time sample
  
  Idx_t=Idx_t+1;        %incremento il contatore della simulazione
  t_curr=t_curr+dt;     %incremento di un dt il tempo della simulazione
  
  %   timeParStr=sprintf( 'Time = %0.5f [s]', t_curr);
  %   hh.text.String=timeParStr;
  %   pause(1e-4)
  %   drawnow
  
  
  if model.compute_recfield==1 & (mod(Idx_t,saveratio)==0) ,
    if(mod(Idx_t,3)==0)
      recsample=m2(sub2ind([Nz Nx],Iz_rec,Ix_rec));
    elseif(mod(Idx_t,3)==1)
      recsample=m0(sub2ind([Nz Nx],Iz_rec,Ix_rec));
    elseif(mod(Idx_t,3)==2)
      recsample=m1(sub2ind([Nz Nx],Iz_rec,Ix_rec));
    end
    
    t_rec=t_rec+1;
    recfield.data(t_rec,:) = recsample(:)';
    recfield.time(t_rec)   = t_curr(:);
    
  end
  
  
end

delete(hh.end)
delete(hh.pause)

if model.compute_recfield==1,
  recfield.data=recfield.data(1:t_rec,:);
  recfield.time=recfield.time(1:t_rec);
end


return
end