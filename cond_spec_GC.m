% Compute the conditional spectral Granger causality for the tri-variate case
% using a non parametric method, spectral matrix factorizzation (Wilson algoritihm)
% parseArgs.m is used to pass additional options
%
% usage:
% [Fxy,Fyx,Fxy_z,Fyx_z,fr,varargout] = cond_spec_GC(x,y,z,varargin)
%
% x, y and z are matrix nxR (n time steps, R realizzations)
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 09 Oct 2009

function [Fxy,Fyx,Fxy_z,Fyx_z,fr,varargout] = cond_spec_GC(x,y,z,varargin)

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
[n,R]=size(x);
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
%inizializzation               
f11=0;f12=0;f13=0;
f21=0;f22=0;f23=0;
f31=0;f32=0;f33=0;

% FT1=zeros(T,n);
% FT2=zeros(T,n);
% FT3=zeros(T,n);

  FT11=zeros(T,n);
  FT12=zeros(T,n);
  FT13=zeros(T,n);
  FT21=zeros(T,n);
  FT22=zeros(T,n);
  FT23=zeros(T,n);
  FT31=zeros(T,n);
  FT32=zeros(T,n);
  FT33=zeros(T,n);

A=nan(3,m*3,R);
C=nan(3,3,R); 

  
for r=1:R
% r/R
% pause(.1)

    if opts.norm>0
    v(:,1)=red2(x(:,r),'anom',anom);
    v(:,2)=red2(y(:,r),'anom',anom);
    v(:,3)=red2(z(:,r),'anom',anom);
    else
    v(:,1)=x(:,r);
    v(:,2)=y(:,r);
    v(:,3)=z(:,r);
    end


%    for i=2:T
% 
%    s0=1/fr(i)/fourier_factor;
%    
%   [FT1(i,:),period(i),scale,coi] = wavelet(v(:,1),dt,PAD,dj,s0,0,mother,param);
%   [FT2(i,:),period(i),scale,coi] = wavelet(v(:,2),dt,PAD,dj,s0,0,mother,param);
%   [FT3(i,:),period(i),scale,coi] = wavelet(v(:,3),dt,PAD,dj,s0,0,mother,param);
% 
%     end
% 
%   f11=f11+mean(FT1.'.*conj(FT1).');
%   f12=f12+mean(FT1.'.*conj(FT2).');
%   f13=f13+mean(FT1.'.*conj(FT3).');
%   f21=f21+mean(FT2.'.*conj(FT1).');
%   f22=f22+mean(FT2.'.*conj(FT2).');
%   f23=f23+mean(FT2.'.*conj(FT3).');
%   f31=f31+mean(FT3.'.*conj(FT1).');
%   f32=f32+mean(FT3.'.*conj(FT2).');
%   f33=f33+mean(FT3.'.*conj(FT3).');
  
  WT = wave(v,'scale',sc(2:end),'dj',dj,'param',param,'dt',dt);
  FT11(2:end,1)=WT(1,1,:);
  FT12(2:end,1)=WT(1,2,:);
  FT13(2:end,1)=WT(1,3,:);
  FT21(2:end,1)=WT(2,1,:);
  FT22(2:end,1)=WT(2,2,:);
  FT23(2:end,1)=WT(2,3,:);
  FT31(2:end,1)=WT(3,1,:);
  FT32(2:end,1)=WT(3,2,:);
  FT33(2:end,1)=WT(3,3,:);
  
  
  f11=f11+FT11;
  f12=f12+FT12;
  f13=f13+FT13;
  f21=f21+FT21;
  f22=f22+FT22;
  f23=f23+FT23;
  f31=f31+FT31;
  f32=f32+FT32;
  f33=f33+FT33;
  
  
  if std(v(:,1))>eps && std(v(:,2))>eps && std(v(:,3))>eps
[~, A(:,:,r), C(:,:,r)]=arfit(v,m,m);

% [w, A1(:,:,r), C1(:,:,r)]=arfit(v(:,1),2,2);
% [w, A2(:,:,r), C2(:,:,r)]=arfit(v(:,2),2,2);
  end
end

%extrapolate spectra at f = 0
%using parametric method for extrap
  A(abs(A)==inf)=nan;
  S = AR_spectrum(nanmean(A,3),nanmean(C,3),fr*dt);S=S*dt;
  use=find(fr>=opts.fextrap & fr>fr(1)); 

  S(1,1,use)=f11(use)/R/df;
  S(1,2,use)=f12(use)/R/df;
  S(1,3,use)=f13(use)/R/df;
  S(2,1,use)=f21(use)/R/df;
  S(2,2,use)=f22(use)/R/df;
  S(2,3,use)=f23(use)/R/df;
  S(3,1,use)=f31(use)/R/df;
  S(3,2,use)=f32(use)/R/df;
  S(3,3,use)=f33(use)/R/df;
  
  
 [SIGMA4,H]=wilson(S,max_iter,tol);
 
 sub=[1 3];
 [SIGMA3,G]=wilson(S(sub,sub,:),max_iter,tol);
 sub=[1 2 3];
 Fyx_z=Granger_multicond(G,H(sub,sub,:),SIGMA3,SIGMA4(sub,sub));

 sub=[2 3];
 [SIGMA3,G]=wilson(S(sub,sub,:),max_iter,tol);
 sub=[2 1 3];
 Fxy_z=Granger_multicond(G,H(sub,sub,:),SIGMA3,SIGMA4(sub,sub));
 
 sub=[1 2];
 [SIGMA2,H2]=wilson(S(sub,sub,:),max_iter,tol);
 [F,Fxy,Fyx] = Granger(SIGMA2,H2,S(sub,sub,:));
% 
% 
%  [SIGMA,H,B1]=wilson(S,max_iter,tol);
%  sub=[1 2 3];
%  Fyx_z=Granger_cond(S(sub,sub,:),H(sub,sub,:),SIGMA(sub,sub));
%  sub=[2 1 3];
%  Fxy_z=Granger_cond(S(sub,sub,:),H(sub,sub,:),SIGMA(sub,sub));

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
plot(fr,[Fxy' Fyx'])
% hold on
% plot(fr,[F' F'+F12'+F21'],'color',[0.4 0.4 0.4]) total and interdependence
% hold on
% plot(fr,[Fp12' Fp21'],'--')  parametric model (of a given order)
legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1]);legend('boxoff')
axis([0 fn 0 max(Fxy+Fyx)])
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
axis([0 fn 0 max(Fxy+Fyx)])
% set(gca,'ylim',yax)
ylabel('spectral G-causality')
xlabel('frequency (time^{-1})')

pause(.1)
end

if opts.test>0
    
h = waitbar(0,'IAAFT test, please wait ...');

for i=1:opts.test
      waitbar(i/opts.test);
%     surrogate=IAAFT(reshape(x,n*R,1));
%     x_s=reshape(surrogate,n,R);
%     Y=repmat(IAAFT(y(:,1)),1,R);
%     [txy(:,i),tyx(:,i),txy_z(:,i),tyx_z(:,i)] =cond_spec_GC(x_s,y,z,'graph','n','disc',opts.disc);
[txy(:,i),tyx(:,i),txy_z(:,i),tyx_z(:,i)] =cond_spec_GC(IAAFT(x),y,z,'graph','n','disc',opts.disc);
    
end
close(h) 
% varargout{1}=G12;
% varargout{2}=G21;


 if opts.test==1
  subplot(subpl(2))
  hold all
  plot(fr,txy,'.-k')
  plot(fr,tyx,'--k')
  legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],...
         [var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1]);legend('boxoff')
  
  subplot(subpl(3))
  hold all
  plot(fr,txy_z,'.-k')
  plot(fr,tyx_z,'--k')
  legend([var1 ' \rightarrow ' var2 ' | ' var3],[var2 ' \rightarrow ' var1 ' | ' var3],...
         [var1 ' \rightarrow ' var2 ' | ' var3],[var2 ' \rightarrow ' var1 ' | ' var3]);legend('boxoff')
  
 else
  subplot(subpl(2))
  hold all
  plot(fr,prctile(txy',95),'.-k')
  plot(fr,prctile(tyx',95),'--k')
  legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],...
         [var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1]);legend('boxoff')
  
  subplot(subpl(3))
  hold all
  plot(fr,prctile(txy_z',95),'.-k')
  plot(fr,prctile(tyx_z',95),'--k')
  legend([var1 ' \rightarrow ' var2 ' | ' var3],[var2 ' \rightarrow ' var1 ' | ' var3],...
         [var1 ' \rightarrow ' var2 ' | ' var3],[var2 ' \rightarrow ' var1 ' | ' var3]);legend('boxoff')
  
 end
end

% 
% function  S = spectral_matrix(s11,s12,s13,s21,s22,s23,s31,s32,s33,use,fr,tr)
% 
% gg=find(fr<1/tr);hh=find(fr>=1/tr);
% a=0;method='cubic';
% 
%   S(1,1,:)=mean(s11(:,use),2);
%   S(1,2,:)=mean(s12(:,use),2);
%   S(1,3,:)=mean(s13(:,use),2);
%   S(2,1,:)=mean(s21(:,use),2);
%   S(2,2,:)=mean(s22(:,use),2);
%   S(2,3,:)=mean(s23(:,use),2);
%   S(3,1,:)=mean(s31(:,use),2);
%   S(3,2,:)=mean(s32(:,use),2);
%   S(3,3,:)=mean(s33(:,use),2);
%   
% 
% for i=1:3
%   for j=1:3
%       
%       if i==j
%       yy(1,:)=S(i,j,hh);
%       y1=fliplr(yy);
%       y2=yy;
%       x1=-fliplr(fr(hh));
%       x2=fr(hh);
%       xi=fr(gg);
%       
%     S(i,j,gg)=interp1([x1 x2],[y1 y2],xi,method);
%     
%       else
%       yy1(1,:)=S(j,i,hh);
%       yy2(1,:)=S(i,j,hh);
%       y1=fliplr(yy1);
%       y2=yy2;
%       x1=-fliplr(fr(hh));
%       x2=fr(hh);
%       xi=fr(gg);
%       
%     S(i,j,gg)=interp1([x1 x2],real([y1 y2]),xi,method)+...
%         sqrt(-1)*interp1([x1 0 x2],imag([y1 0 y2]),xi,method);
%     
%       end
% %        a=a+1;
% %   subplot(3,3,a)
% %   ff(:,1)=S(i,j,:);semilogx(fr,[real(ff) imag(ff)],'.-')     
%   end
%   
% end
% 
% end
% end
