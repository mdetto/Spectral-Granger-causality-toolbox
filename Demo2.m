%this Demo2 show three examples of bivariate interactions
%%
clear
close all
addpath('./ARFIT/');

%options
n=512;
R=50;
dlag=5;
disc=50;
test=0;
param=15;
alfa=1;


w=[0 0];      % zero mean processes
C=[1 0.0; 0.0 1]; % covariance matrix (canonical form)

A1=[0.55 0.0;0.0 0.55];  % 1-st order coefficients of the AR model
A2=[-0.8 0.0;0.0 -0.8];  % 2-st order coefficients of the AR model

v = arsim(w,[A1 A2],C,n*R,n);
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),...
    'disc',disc,'subpl',331,'test',test,'param',param,...
     'method','wave');

[Y lag]=xcov(v(:,1),v(:,2),dlag,'coeff');
subplot(337)
plot(lag,Y)
ylabel('cross-correlation')
xlabel('lag')
pause(.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=[0.55 0.2;0.0 0.55];
A2=[-0.8 0.0;0.0 -0.8];

v = arsim(w,[A1 A2],C,n*R,n);
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),...
    'disc',disc,'subpl',332,'test',test,'param',param);

[Y lag]=xcov(v(:,1),v(:,2),dlag,'coeff');
subplot(338)
plot(lag,Y)
xlabel('lag')
pause(.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%

A1=[0.55 0.2;0.1 0.55];
A2=[-0.8 0.0;0.0 -0.8];

v = arsim(w,[A1 A2],C,n*R,n);
spec_GC(reshape(v(:,1),n,R),reshape(v(:,2),n,R),...
    'disc',disc,'subpl',333,'test',test,'param',param);

[Y lag]=xcov(v(:,1),v(:,2),dlag,'coeff');
subplot(339)
plot(lag,Y)
xlabel('lag')
pause(.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%set titles and axis of spectral-G
subplot(331)
title(' no interaction')
subplot(332)
title(' simple interaction')
subplot(333)
title(['interaction and feedback'])

subplot(336)
h=get(gca,'ylim');
subplot(334)
set(gca,'ylim',h)
subplot(335)
set(gca,'ylim',h)

subplot(338)
h=get(gca,'ylim');
subplot(337)
set(gca,'ylim',h)
subplot(339)
set(gca,'ylim',h)


