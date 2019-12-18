% This demo shows the use of the bivariate spectral Granger Causality
% using 1)parametric method, 2) Fourier transform 3) Wavelet
% Two systems of 2-nd order with interaction between X and Y are simulated
%
%  x(i)= 0.55x(i-1) - 0.8x(i-2) + 0.2y(i-1) 
%  y(i)= 0.55y(i-1) - 0.8y(i-2) 
%%
clear
close all
addpath('./ARFIT/');

n=100;  % number of points 
R=500;  % number of realizzations
test= 0; % number od montecarlo simulation for significance usinf IAAFT data
A1=[0.55 0.4;0.0 0.55];  % 1-st order coefficients of the AR model
A2=[-0.8 0.0;0.0 -0.8];  % 2-st order coefficients of the AR model

w=[0 0];
C=[0.7 0.0; 0.0 0.7]; % covariance matrix (canonical form)


v = arsim(w,[A1 A2],C,n*R,n);

spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),'disc',50,'subp',331,'method','para');pause(.1)
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),'disc',50,'subp',332,'method','four');pause(.1)
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),'disc',50,'subp',333,'method','wave');

%set titles and axis of spectral-G
subplot(331)
title('parametric')
subplot(332)
title('Fourier')
subplot(333)
title('Wavelet')



