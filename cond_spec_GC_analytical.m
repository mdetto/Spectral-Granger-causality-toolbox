% Compute the conditional spectral Granger causality analytically
% passing the parametyer matrix A of a AR model and covariance matrix C
% parseArgs.m is used to pass additional options
%
% usage:
% [Fxy_z,Fyx_z,fr] = cond_spec_GC_analytical(A,C,varargin)
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 09 Oct 2009

function [Fxy_z,Fyx_z,fr] = cond_spec_GC_analytical(A,C,varargin)

%default options
opts=struct('mother','Morlet',... %mother wavelet (need to be complex)
            'param',15,...        %wavelet base parameter
            'dt',1,...            %temporal step
            'max_iter',20, ...    %max iterations for Wilson algorithim
            'tol',1e-6, ...       %tolerance for Wilson algorithim
            'fextrap',0, ...      %limit for spectra extrapolation 
            'order',2,...         %order of the AR model for arfit.m 
            'disc',100,...        %number of intervals of the Fourier domain 
            'norm',0,...          %normilize each trial 
            'anom',0,...          %seasonal de-trend, anom = number of season
            'test',0,...          %number of IAAFT data for sign-test
            'graph','y',...       %plot results y/n
            'subpl',321,...       %spectral will be plotted in 321, GC in 323, Cond-GC in 325
            'vname',0);           %pass the name of the varaibles as cell array    
opts=parseArgs(varargin,opts);


%constant and parameters
m=opts.order;
dt=opts.dt;
fn=1/2/dt;             %Niquist frequency
dth=pi/(opts.disc-1);  %interval of discretizzation of angular frequency
yax=[0 5];

max_iter=opts.max_iter;
tol=opts.tol;
anom=opts.anom;

theta=0:dth:pi;% discretizzation of fourier space (in angular frequency)
T=length(theta);
S=zeros(2,2,T);
fr=theta/2/pi/dt;

%paramentrs and costant wavelet
dj=0.1;      
Cdelta = 0.776;
psi0=pi^-0.25;
PAD=1;
J=-1;
mother=opts.mother;
param=opts.param;
fourier_factor = (4*pi)/(param + sqrt(2 + param^2));
df=Cdelta*2*pi*log(2)/dt;   % spectral interval 
sc=1./fr/fourier_factor;


%using parametric method for extrap

  [S,H,SIGMA] = AR_spectrum(A,C,fr*dt);S=S*dt;
 
 sub=[1 2 3];
 Fyx_z=Granger_cond(S(sub,sub,:),H(sub,sub,:),SIGMA(sub,sub));
 
 sub=[2 1 3];
 Fxy_z=Granger_cond(S(sub,sub,:),H(sub,sub,:),SIGMA(sub,sub));
 

if opts.graph == 'y'

    %variables names
    if ~iscell(opts.vname)
        var1='x';
        var2='y';
        var3='z';
    else
        var1=opts.vname{1};
        var2=opts.vname{1};
        var3=opts.vname{1};
    end

    
    subpl(1)=opts.subpl;
    ncol=num2str(subpl);
    subpl(2)=subpl(1)+str2num(ncol(2));
    subpl(3)=subpl(2)+str2num(ncol(2));
subplot(subpl(1))
s11(:,1)=S(1,1,:);
s22(:,1)=S(2,2,:);
s33(:,1)=S(3,3,:);
plot(fr,[s11 s22 s33])
legend(var1,var2,var3);legend('boxoff')
axis([0 fn 0 max(s11+s22)])
ylabel('wavelet spectra')
xlabel('frequency (time^{-1})')
subplot(subpl(2))
plot(fr,fr*0)
% hold on
% plot(fr,[F' F'+F12'+F21'],'color',[0.4 0.4 0.4]) total and interdependence
% hold on
% plot(fr,[Fp12' Fp21'],'--')  parametric model (of a given order)
legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1]);legend('boxoff')
axis([0 fn 0 1])
% set(gca,'ylim',yax)
ylabel('spectral G-causality')
xlabel('frequency (time^{-1})')

subplot(subpl(3))
plot(fr,[Fxy_z' Fyx_z'])
% hold on
% plot(fr,[F' F'+F12'+F21'],'color',[0.4 0.4 0.4]) total and interdependence
% hold on
% plot(fr,[Fp12' Fp21'],'--')  parametric model (of a given order)
legend([var1 ' \rightarrow ' var2 ' | ' var3],[var2 ' \rightarrow ' var1 ' | ' var3]);legend('boxoff')
axis([0 fn 0 1.2*max(Fxy_z+Fyx_z)])
% set(gca,'ylim',yax)
ylabel('spectral G-causality')
xlabel('frequency (time^{-1})')

pause(.1)
end

