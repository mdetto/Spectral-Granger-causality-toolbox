%wavelet.m:   Computes wavelet spectral matrix Morlet Continous Wavelet Transform

%
%       [WT,scale] = wavelet(s);
%
%       [WT,scale] = wavelet(s,options)

%
%   		s       = n x m time series  
%           WT      =  spectral density matrix [S11 S12; S21 S22];
%           scale   = frequency (2*pi*n/dt) - optional          
%
%   options:
%
%   'dj' = the spacing between discrete scales. Default is 0.25.
%          A smaller # will give better scale resolution, but be slower to plot.
%
%   's0' = the smallest scale of the wavelet.  Default is 2*DT.
%
%   'J1' = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%          to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%   'mother' = the mother wavelet function. The choices are 'MORLET',
%              'PAUL', or 'DOG' - default 'Morlet'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%
%      'dt'      =  interval time of the signal - default 1
%      'graph'    =  'y' 'n' - default 'N'
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 5 April 2010
% -------------------------------------------------------------------------
%   Copyright (C) 2007-2010, Matteo Detto
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


function [WT,scale] = wave(v,varargin)


%default options
opts=struct('mother','Morlet',...  
            's0',2, ...    
            'PAD',0, ...    
            'param',[6 1], ...    
            'angle',0,...
            'J1',-1, ...    
            'dj',0.15, ... 
            'dt',1,...
            'scale',0,...
            'graph','n');
opts=parseArgs(varargin,opts);

[n,m]=size(v);
s0=opts.s0;
dj=opts.dj;

if opts.J1==-1 && opts.scale(1)==0
    J1=floor((log(n/2/s0)/log(2))/dj);
    scale(:,1) = s0*2.^((0:J1)*dj);
elseif opts.J1>0 && opts.scale(1)==0
    J1=fix((log(opts.J1/s0/opts.dx)/log(2))/dj);
    scale(:,1) = s0*2.^((0:J1)*dj);
elseif opts.scale(1)>0
    scale(:,1)=opts.scale/opts.dt;
end

S=length(scale);
npuls  = floor((n-1)/2);
k(:,1)  = 2*pi/n*[ 0:npuls  (npuls-n+1):-1 ];
K=length(k);

k0=opts.param(1);

gap=find(isnan(v));
if isempty(gap)
FT=zeros(K,m);
for i=1:m
 FT(:,i)=fft(v(:,i));
end

FS=zeros(m,m,K);
for i=1:m
      for j=1:m
         FS(i,j,:)=FT(:,i).*conj(FT(:,j));
      end
end

%%%%%%%%%%%%%%
else
  FS=zeros(m,m,K);
B=nan(n,1);
for i=1:m
      for j=1:m
  %%%%%%%%%%%%        
  if j==i 
   for r=1:n/2+1
    B(r)=nanmean(v(1:n-r+1,i).*v(r:n,j))*n;
   end
    B(n/2+2:n)=B(n/2:-1:2);
    gap=find(isnan(B));
    if ~isempty(gap)
    use=find(~isnan(B));
    B(gap)=interp1(use,B(use),gap);
    end
    FS(i,j,:)=ifft(B);
    %%%%%%%%%%%
  elseif j>i

    for r=1:n/2
        B(r)=nanmean(v(1:n-r+1,i).*v(r:n,j))*n;
        B(n-r+1)=nanmean(v(1:n-r,j).*v(r+1:n,i))*n;
    end
    gap=find(isnan(B));
    if ~isempty(gap)
    use=find(~isnan(B));
    B(gap)=interp1(use,B(use),gap);
    end
    FS(i,j,:)=ifft(B);
%%%%%%%%%%%%%
  elseif j<i
      FS(i,j,:)=conj(FS(j,i,:));
  end


      end
end
end
WT=zeros(m,m,S);
for s=1:S

%     H=abs(psi(scale(s)*k,k0)).^2;
    H=exp(-(scale(s)*k-k0).^2);
    %numerical normalizzation
    norm=sum(H)*n;
    %analitical normalizzation
    %norm=1/sqrt(pi)*n;
    
    for i=1:m
      for j=1:m
         W(:,1)=FS(i,j,:);
         WT(i,j,s)=sum(W.*H)/norm;
      end
    end
end
    
% S11(:,1)=WT(1,1,:);
% loglog(scale,S11);
% pause
%convert wavelet scale in Fourier period
fourier_factor=4*pi/(k0+sqrt(2+k0^2));
scale=scale*opts.dt*fourier_factor;


