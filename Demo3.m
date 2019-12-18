% This demo shows the use of the tri-variate conditional spectral Granger Causality
%
% Two systems of 2-nd order with/without interaction between X, Y and Z are simulated
%
%  x(i)= 0.55x(i-1) - 0.8x(i-2) + by(i-1) 
%  y(i)= 0.55y(i-1) - 0.8y(i-2) 
%  z(i)= 0.55z(i-1) - 0.8y(i-2) 
%
% where b = 0.2 for system 1 and b = 0.0  for system 2
%%
clear
close
addpath('./ARFIT/');

%---------------------------------------
%options
disc=50;
param=15;
norm=1;

R=100;  % number of realizzations
test= 0; % number od montecarlo simulation for significance usinf IAAFT data

%   y > x | z;  
%paramentrs and costant
n=200;                   %# of points
A1=[0.0 0.0 0.4
    0.0 0.5 0.4
    0.0 0.0 0.55];        %AR1 coeff matrx
A2=[0.0  0.4  0.0
    0.0  -0.8  0.0
    0.0  0.0 -0.8];       %AR2 coeff matrx

% A1=[0.55 0.0 0.0
%     0.0 0.55 0.0
%     0.0 0.0 0.55];        %AR1 coeff matrx
% A2=[-0.8  0.4  0.0
%     0.0  -0.8  0.0
%     0.0  0.0 -0.8];       %AR2 coeff matrx


w=zeros(3,1);

C=eye(3);  %canonical form
%add some cross-correlation to test the normalization
% C(2,1)=0.25;C(1,2)=0.25;  

v = arsim(w,[A1 A2],C,n*R,n);
x=v(:,1);
y=v(:,2);
z=v(:,3);


cond_spec_GC(reshape(x,n,R),reshape(y,n,R),reshape(z,n,R),...
    'disc',disc,'subpl',321,'norm',norm);


A2(1,2)=0;
v = arsim(w,[A1 A2],C,n*R,n);
x=v(:,1);
y=v(:,2);
z=v(:,3);

cond_spec_GC(reshape(x,n,R),reshape(y,n,R),reshape(z,n,R),...
    'disc',disc,'subpl',322,'norm',norm);

%set titles
subplot(321)
title('y has direct interaction on x')
subplot(322)
title(['y has not direct interaction';...
       'on x (but mediated by z)    '])
