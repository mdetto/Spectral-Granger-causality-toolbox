%Compute spectral Granger causality analytically for a give AR model  
%with coeff matrix A and noise matrix C
% parseArgs.m is used to pass additional options
%
% [F12,F21,fr,varargout] = spec_GC(x,y,varargin);
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 09 Oct 2009

function [F12,F21,F,fr,varargout] = spec_GC_analytical(A,C,varargin)

%default options
opts=struct('dt',1,...            %temporal step
            'disc',100,...        %number of intervals of the Fourier domain 
            'graph','y',...       %plot results y/n
            'subpl',321,...       %subplot results y/n
            'vname',['x' 'y']);   %variable names
opts=parseArgs(varargin,opts);


%constant and parameters
dt=opts.dt;

fn=1/2/dt;             %Niquist frequency
dth=pi/(opts.disc-1);  %interval of discretizzation of angular frequency
yax=[0 5];

theta=0:dth:pi;% discretizzation of fourier space (in angular frequency)
T=length(theta);
fr=theta/2/pi/dt;

%using parametric method
  [S,H,SIGMA] = AR_spectrum(A,C,fr*dt);
% s11(:,1)=S(1,1,:);
% plot(fr,s11,'.-')
% pause

%   S=S*dt;
%   H=H*sqrt(dt);

 [F,F12,F21] = Granger(SIGMA,H,S);F21

if opts.graph == 'y'

    subpl(1)=opts.subpl;
    ncol=num2str(subpl);
    subpl(2)=subpl(1)+str2num(ncol(2));
    if iscell(opts.vname)
    var1=opts.vname{1};
    var2=opts.vname{2};
    else
    var1=opts.vname(1);
    var2=opts.vname(2);
    end
    
subplot(subpl(1));

s11(:,1)=S(1,1,:);
s22(:,1)=S(2,2,:);
s12(:,1)=S(1,2,:);
coh=abs(s12).^2./s11./s22;

plot(fr,[s11 s22])
axis([0 fn 0 max(s11+s22)])
ylabel('wavelet spectra')
legend(var1,var2);legend('boxoff')
subplot(subpl(2));
plot(fr,[F12' F21'])
hold on
plot(fr,[F' F'+F12'+F21'],'color',[0.4 0.4 0.4]) %total and interdependence
% hold on
% plot(fr,[Fp12' Fp21'],'--')  parametric model (of a given order)
legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],'inter','total');legend('boxoff')
axis([0 fn 0 max(F12+F21)])
% set(gca,'ylim',yax)
ylabel('spectral G-causality')
xlabel('frequency (time^{-1})')
%  finfo=dir('*.fig');
%  fnum = num2str(length(finfo)+1);
%  saveas(gcf,['spec_GC_n' fnum '.fig'])
pause(.1)
end



