function [fact]=seisplot2(varargin)
% [fact]=seisplot2(datain,t,tr,scal,pltflg,scfact,colour,clip)
%
% function for plotting seismic traces
%
% INPUT
% datain  - input matrix of seismic traces
% t       - time axis
% tr      - trace axis
% scal    - 1 for global max, 0 for global ave, 2 for trace max
% pltflg  - 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
%           2 plots wiggle traces only, 3 imagesc gray, 4 pcolor gray
% scfact  - scaling factor
% colour  - trace colour, default is black
% clip    - clipping of amplitudes (if <1); default no clipping
%
% OUPTPUT
% fact    - factor that matrix was scaled by for plotting
% if you want to plot several matrices using the same scaling factor,
% capture 'fact' and pass it as input variable 'scal' to future matrices
% with 'scfact' set to 1
%
% DSI customized VSP processing software
% written by G. Perron January, 1998

%Revision 3.2  2008/12/15 gianca
%added clip  
%
%Revision 3.1  2006/03/30 gianca
%added varagin and other controls  
%
%$Id: seisplot.m,v 3.0 2000/06/13 19:19:15 gilles Exp $
%$Log: seisplot.m,v $
%Revision 3.0  2000/06/13 19:19:15  gilles
%Release 3
%
%Revision 2.0  1999/05/21 18:44:00  mah
%Release 2
%
%Revision 1.5  1999/05/21 17:56:09  mah
%checked in file which wasn't checked in properly
%
%Revision 1.4  1999/05/21 17:36:41  mah
%checking in file which wasn't checked in properly
%
%Revision 1.3  1999/01/14 21:08:18  perron
%Fixing a typo in the previous version
%
%Revision 1.2  1999/01/14 20:43:10  perron
%Adding tstart variable to take into account subset variables when computing sample range to plot
%
%Revision 1.1  1999/01/06 19:08:35  kay
%Initial revision
%
%Copyright (C) 1998 Seismology and Electromagnetic Section/
%Continental Geosciences Division/Geological Survey of Canada
%
%This library is free software; you can redistribute it and/or
%modify it under the terms of the GNU Library General Public
%License as published by the Free Software Foundation; either
%version 2 of the License, or (at your option) any later version.
%
%This library is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%Library General Public License for more details.
%
%You should have received a copy of the GNU Library General Public
%License along with this library; if not, write to the
%Free Software Foundation, Inc., 59 Temple Place - Suite 330,
%Boston, MA  02111-1307, USA.
%
%DSI Consortium
%Continental Geosciences Division
%Geological Survey of Canada
%615 Booth St.
%Ottawa, Ontario
%K1A 0E9
%
%email: dsi@cg.nrcan.gc.ca

if nargin==0; help seisplot2; return; end

a=(varargin{1});
[M,N]=size(a);

if M<N, 
    disp('Warning: traces # > time sample #'); 
%    ok=input('Continue anyway (y/n)? -> ','s');
%    if ok~='y', return, end
end

% defaults
t  = (0:M-1); 
traces = 1:N; 
scal   = 1;
pltflg = 0;
scfact = 1;
colour = 'k';
clip   =  1;

lv=length(varargin);
for k=lv+1:8; varargin{k}=[];end

if ~isempty(varargin{2}),  t=varargin{2};       end;
if ~isempty(varargin{3}),  traces=varargin{3};  end;
if ~isempty(varargin{4}),  scal=varargin{4};    end;
if ~isempty(varargin{5}),  pltflg=varargin{5};  end;
if ~isempty(varargin{6}),  scfact=varargin{6};  end;
if ~isempty(varargin{7}),  colour=varargin{7};  end;
if ~isempty(varargin{8}),  clip=varargin{8};  end;

Na1=size(a,1);
Na2=size(a,2);

if clip<1;
    for k1=1:N;
        aa=a(:,k1)-mean(a(:,k1));
        [abin,aval]=hist(abs(aa),200);
        aas=cumsum(abin)/sum(abin);
        II=min(find(aas>clip));
        aa(find(abs(aa)>aval(II)))=sign(aa(find(abs(aa)>aval(II)))).*aval(II);
        a(:,k1)=aa+mean(a(:,k1));
    end
end

if scal==1
   fact=max(max(abs(a)));
elseif scal==0
   fact=max(mean(abs(a)));
elseif scal==2
   fact=1;
   fact1=max(abs(a));
   fact1(find(fact1==0))=1;
   a=a./repmat(fact1,Na1,1);
else
   fact=scal;
end
fact=fact./scfact;
a=a./fact;

if pltflg==3
  imagesc(traces,t,a);
  colormap(gray)
  colorbar
  xlabel('Trace Axis')
  ylabel('Time (s)')
  return;
end %if

if pltflg==4
  pcolor(traces,t,a); shading flat
  colormap(gray)
  colorbar
  xlabel('Trace Axis')
  ylabel('Time (s)')
  set(gca,'ydir','reverse')
  return;
end %if

b=find(a<0);
c=a;
c(b)=0;
inc=mean(abs(gradient(traces)));
%[xmat,ymat]=meshgrid(traces,t1:smp:t2);
%c=c.*inc+xmat;
%c(1,:)=xmat(1,:);
%c(round(t2/smp+1-t1/smp),:)=xmat(1,:);
%a=a.*inc+xmat;

[xmat,ymat]=meshgrid(traces,t);
c=c.*inc+xmat;
c(1,:)=xmat(1,:);
c(end,:)=xmat(1,:);
a=a.*inc+xmat;

if pltflg==0
   h=fill(c,ymat,colour);
   set(h,'edgecolor','none')
   hold on
   plot(a,ymat,colour);
   hold off
elseif pltflg==1
   h=fill(c,ymat,colour);
   set(h,'edgecolor','none')
elseif pltflg==2
   plot(a,ymat,colour);
end

set(gca,'ydir','reverse')
%axis([min(traces) max(traces) min(t) max(t)])
axis([min(min(a)) max(max(a)) min(t) max(t)])

xlabel('Trace axis')
ylabel('Time (s)')
set(gca,'ygrid','on')
