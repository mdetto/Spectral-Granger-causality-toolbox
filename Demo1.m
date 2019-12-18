% This demo shows the use of the bi-variate spectral Granger Causality
% using a non parametric method. 
% Two systems of 2-nd order with/without interaction between X and Y are simulated
%
%  x(i)= 0.55x(i-1) - 0.8x(i-2) + by(i-1) + e1(i)
%  y(i)= 0.55y(i-1) - 0.8y(i-2) + e2(i)
%
% where b = 0.4 for system 1 and b = 0.0  for system 2, e1 and e2 are
% correlated
%%
clear
close all
addpath('./ARFIT/');

% case 1 b = 0.4
b1= 0.4;
% case 2 b = 0.0 (no interactions)
b2= 0.0;

n=25;  % number of points 
R=25;  % number of realizzations
test= 0; % number od montecarlo simulation for significance usinf IAAFT data
disc=round(n/2);
w=[0 0];
C=[0.7 0.5; 0.5 0.7]; % covariance matrix


% simulate process 1
A1=[0.55 b1; 0.0 0.55];  % 1-st order coefficients of the AR model
A2=[-0.8 0.0;0.0 -0.8];  % 2-st order coefficients of the AR model
v = arsim(w,[A1 A2],C,n*R,n);

plot(v(1:100,:))
legend('x','y')
title('who causes who ?')
pause
close

%compte GC and plot
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),...
    'disc',disc,'subp',321,'test',test,'param',25,'method','wave','fextrap',0.02);


% simulate process 2
A1=[0.55 b2; 0.0 0.55];  % 1-st order coefficients of the AR model
A2=[-0.8 0.0;0.0 -0.8];  % 2-st order coefficients of the AR model
v = arsim(w,[A1 A2],C,n*R,n);

%compute GC and plot
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),'disc',disc,'subp',322,'test',test,'fextrap',0.02);

%set titles and axis of spectral-G
subplot(321)
title([' b = ' num2str(b1)])
subplot(322)
title([' b = ' num2str(b2)])
subplot(323)
h=get(gca,'ylim');
subplot(324)
set(gca,'ylim',h)


