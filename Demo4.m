% This demo shows the use of the multi-variate conditional spectral Granger Causality 
% look the AR coefficient matrix to understand direct interactions
%%
clear
close all
addpath('./ARFIT/');

%options
n=100;  % number of points 
R=100;  % number of realizzations
disc=round(n/2);
test=0;
param=15;
fextrap=0;
norm=1;



%   y > x | z;  
%coefficient matrix of the AR model
A1=[0.8 0.2 0.0
    0.0 0.9 0.0
    0.0 0.0 0.5];        %AR1 coeff matrx
A2=[-0.5  0.0  0.0
    0.0  -0.8  0.0
    0.0  0.0 -0.2];       %AR2 coeff matrx

% A1=[0.55 0.0 0.0
%     0.0 0.55 0.0
%     0.0 0.0 0.55];        %AR1 coeff matrx
% A2=[-0.8  0.4  0.0
%     0.0  -0.8  0.0
%     0.0  0.0 -0.8];       %AR2 coeff matrx


w=zeros(3,1);

% C=eye(3);  %canonical form
% %add some cross-correlation to test the normalization
% C(2,1)=0.25;C(1,2)=0.25;  
C=[0.3 0 0
   0   1 0  
   0   0 0.2];

%simulate AR
v = arsim(w,[A1 A2],C,n*R,n);


%the variables parsed into the functions as cell array
X{1}=reshape(v(:,1),n,R);
X{2}=reshape(v(:,2),n,R);
X{3}=reshape(v(:,3),n,R);
% X{4}=randn(n,R); %dummies variables independent from x, y and z
% X{5}=randn(n,R); %dummies variables independent from x, y and z
% X{6}=randn(n,R); %dummies variables independent from x, y and z
vname{1}='X';
vname{2}='Y';
vname{3}='Z';

multi_spec_GC(X,'fextrap',fextrap,'disc',disc,'norm',norm,'vname',vname);



