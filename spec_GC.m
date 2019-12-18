%Compute spectral Granger causality for the bivariate 
%using a non parametric method, spectral matrix factorizzation (Wilson algoritihm)
% parseArgs.m is used to pass additional options
%
% [F12,F21,fr,varargout] = spec_GC(x,y,varargin);
%
% x and y are matrix nxR, n time steps, R realizzations
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 09 Oct 2009
% 
%default options
% opts=struct('method','wave',...   % or 'four' or 'para'
%             'A',0,...             % if the AR model is known (for method para)
%             'SIGMA',0,...         % if error covariance matrix is know (for method para)  
%             'mother','Morlet',... %mother wavelet (need to be complex)
%             'param',15,...        %wavelet base parameter
%             'dt',1,...            %temporal step
%             'max_iter',20, ...    %max iterations for Wilson algorithim
%             'tol',1e-6, ...       %tolerance for Wilson algorithim
%             'fextrap',0, ...      %limit for spectra extrapolation 
%             'order',2,...         %order of the AR model for arfit.m 
%             'disc',100,...        %number of intervals of the Fourier domain 
%             'norm',0,...          %normilize each trial
%             'anom',[0 0],...      %remove periodicity anom
%             'test',0,...          %number of IAAFT data for sign-test
%             'graph','y',...       %plot results y/n
%             'subpl',321,...       %subplot results y/n
%             'vname',['x' 'y']);   %variable names
        
function [F12,F21,F,fr,varargout] = spec_GC(x,y,varargin)

%default options
opts=struct('method','wave',...   % or 'four'
            'A',0,...             % if the AR model is known (for method para)
            'SIGMA',0,...         % if error covariance matrix is know (for method para)  
            'mother','Morlet',... %mother wavelet (need to be complex)
            'param',15,...        %wavelet base parameter
            'dt',1,...            %temporal step
            'max_iter',20, ...    %max iterations for Wilson algorithim
            'tol',1e-6, ...       %tolerance for Wilson algorithim
            'fextrap',0, ...      %limit for spectra extrapolation 
            'order',2,...         %order of the AR model for arfit.m 
            'disc',100,...        %number of intervals of the Fourier domain 
            'norm',0,...          %normilize each trial
            'anom',[0 0],...      %remove periodicity anom
            'test',0,...          %number of IAAFT data for sign-test
            'graph','y',...       %plot results y/n
            'subpl',321,...       %subplot results y/n
            'vname',['x' 'y']);   %variable names
opts=parseArgs(varargin,opts);


%constant and parameters
m=opts.order;
dt=opts.dt;
[n,R]=size(x);
[n2,R2]=size(y);
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
f11=0;f12=0;
f21=0;f22=0;

 if strcmp(opts.method,'wave')
   FT11=zeros(T,1);
   FT21=zeros(T,1);
   FT12=zeros(T,1);
   FT22=zeros(T,1);
%    FT1=zeros(T,n);
%    FT2=zeros(T,n);

 elseif strcmp(opts.method,'four')
   FT1=zeros(T,1);
   FT2=zeros(T,1);
 elseif strcmp(opts.method,'para')
   A =nan(2,m*2,R);C =nan(2,2,R); 
   A1=nan(m,R);C1=nan(R,1); 
   A2=nan(m,R);C2=nan(R,1); 
 end
  
for r=1:R
 R-r

    if R2==1  %if y has only one relaizzation
        r2=1;
    else
        r2=r;
    end
    
    if opts.norm>0
    v(:,1)=red2(x(:,r),'anom',anom(1),'pad','y');
    v(:,2)=red2(y(:,r2),'anom',anom(2),'pad','y');
    else
    v(:,1)=x(:,r);
    v(:,2)=y(:,r2);
    end

if opts.method=='wave'
    
%    for i=2:T
%     s0=1/fr(i)/fourier_factor;
%     [FT1(i,:),period(i),scale,coi] = wavelet(v(:,1),dt,PAD,dj,s0,0,mother,param);
%     [FT2(i,:),period(i),scale,coi] = wavelet(v(:,2),dt,PAD,dj,s0,0,mother,param);
%    end
% 
%   f11=f11+mean(FT1.*conj(FT1),2);
%   f12=f12+mean(FT1.*conj(FT2),2);
%   f21=f21+mean(FT2.*conj(FT1),2);
%   f22=f22+mean(FT2.*conj(FT2),2);

  WT = wave(v,'scale',sc(2:end),'dj',dj,'param',param,'dt',dt);
  FT11(2:end,1)=WT(1,1,:);
  FT21(2:end,1)=WT(1,2,:);
  FT12(2:end,1)=WT(2,1,:);
  FT22(2:end,1)=WT(2,2,:);
  
  f11=f11+FT11;
  f12=f12+FT21;
  f21=f21+FT12;
  f22=f22+FT22;
  
elseif  opts.method=='four'
    
     for i=2:T
        period(i)=1/fr(i);
        FT1(i,1) = fourier(v(:,1),theta(i));
        FT2(i,1) = fourier(v(:,2),theta(i));
     end

  f11=f11+FT1.*conj(FT1)/n;
  f12=f12+FT1.*conj(FT2)/n;
  f21=f21+FT2.*conj(FT1)/n;
  f22=f22+FT2.*conj(FT2)/n;
  coi=1/fr(2);
end

%parametric fit (predefined m-order)
 [w, A(:,:,r), C(:,:,r)]=arfit(v,m,m);
 [w, A1(:,r), C1(r)]=arfit(v(:,1),m,m);
 [w, A2(:,r), C2(r)]=arfit(v(:,2),m,m);

end

%extrapolate spectra at f = 0
%using parametric method
  A(abs(A)==inf)=nan;
  [S,H,SIGMA] = AR_spectrum(nanmean(A,3),nanmean(C,3),fr*dt);
  S=S*dt;
  H=H*sqrt(dt);
%   [SIGMA,H,B1,B2]=wilson(S,max_iter,tol); %spectral matrix factorizzation
%   [F,F12,F21] = Granger(SIGMA_p,H,S);
%check here how good is the model fitting
% figure
% subplot(321)
% s11(:,1)=S(1,1,:);
% loglog(fr,s11);hold all
% loglog(fr,f11/R/df)
% subplot(325)
% s22(:,1)=S(2,2,:);
% loglog(fr,s22);hold all
% loglog(fr,f22/R/df)
% pause

if opts.method=='wave' | opts.method=='four'
%extrapolate
use=find(fr>=opts.fextrap & fr>fr(1)); 

  S(1,1,use)=f11(use)/R/df;
  S(1,2,use)=f12(use)/R/df;
  S(2,1,use)=f21(use)/R/df;
  S(2,2,use)=f22(use)/R/df;
  
  display('running Wilson algorithm, please wait...')
  [SIGMA,H,B1,B2]=wilson(S,max_iter,tol); %spectral matrix factorizzation

  end
      
   
  [F,F12,F21] = Granger(SIGMA,H,S); %spectral G-causality
   

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
plot(fr,[F12' F21']);hold all %directional
% plot(fr,F,'--','color',[0.4 0.4 0.4]); %interdependence
plot(fr,[F+F12+F21],'color',[0.4 0.4 0.4])% total
% hold on
% plot(fr,[Fp12' Fp21'],'--')  parametric model (of a given order)
legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1]);legend('boxoff')
%legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],'inter','total');legend('boxoff')
axis([0 fn 0 max(F12+F21)])
% set(gca,'ylim',yax)
ylabel('spectral G-causality')
xlabel('frequency (time^{-1})')
%  finfo=dir('*.fig');
%  fnum = num2str(length(finfo)+1);
%  saveas(gcf,['spec_GC_n' fnum '.fig'])
pause(.1)
end

if opts.test>0
    
    h = waitbar(0,'IAAFT test, please wait...');
for i=1:opts.test
      waitbar(i/opts.test);
%     surrogate=IAAFT(reshape(x,n*R,1));
%     X=reshape(surrogate,n,R);
%     surrogate=IAAFT(reshape(y,n*R,1));
%     Y=reshape(surrogate,n,R);
  if n>=50
    X=IAAFT(x);
    Y=IAAFT(y);
  else %series too short for IAAFT
    X=x(randperm(n),:);
    Y=y(randperm(n),:);
  end

% C = chol(SIGMA);                    % R is upper triangular
%  for d=1:R
%   vr = [x(:,d) y(:,d)]*inv(C);
%   ir(:,1)=IAAFT(vr(:,1));
%   ir(:,2)=IAAFT(vr(:,2));
%   hr=ir*C;
%   X(:,d)=hr(:,1);
%   Y(:,d)=hr(:,2);
%  end
   [G12(:,i),G21(:,i),G(:,i)] =spec_GC(X,Y,'graph','n','disc',opts.disc,...
       'method',opts.method,'order',opts.order,'param',opts.param);
    
end
   close(h) 
varargout{1}=G12;
varargout{2}=G21;

if opts.test==1
   hold all
   plot(fr,G21,'.-k')
   plot(fr,G12,'-k')
   legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],...
        't21(.95)','t12(.95)');legend('boxoff')
else
    hold all
    plot(fr,prctile(G21',95),'.-k')
    plot(fr,prctile(G12',95),'-k')
%     plot(fr,prctile(G',95),'--k')
%     legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],...
%         't21(.95)','t12(.95)','t_{inter}');legend('boxoff')
    legend([var1 ' \rightarrow ' var2],[var2 ' \rightarrow ' var1],...
        't21(.95)','t12(.95)');legend('boxoff')

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
