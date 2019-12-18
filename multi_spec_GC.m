% Compute spectral Granger causality for the multi-variate 
% using a non parametric method, spectral matrix factorizzation (Wilson algoritihm)
% parseArgs.m is used to pass additional options
%
% [Fyx,Fxy,Fyx_z,varargout] = multi_spec_GC(x,varargin)
% the variables parsed into the functions as cell array
% X{1}=reshape(x,n,R);
% X{2}=reshape(y,n,R);
% X{3}=reshape(z,n,R);
% x, y and z are vectors reshaped into matrix nxR, n time steps, R realizzations
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 1 April 2011

function [Fyx,Fxy,Fyx_z,fr] = multi_spec_GC(x,varargin)

%default options
opts=struct('mother','Morlet',... %mother wavelet (need to be complex)
            'param',15,...        %wavelet base parameter
            'dt',1,...            %temporal step
            'max_iter',30, ...    %max iterations for Wilson algorithim
            'tol',1e-6, ...       %tolerance for Wilson algorithim
            'fextrap',0, ...      %limit for spectra extrapolation 
            'extrap_method','param', ...      %method of extrapolation (param or interp) 
            'order',2,...         %order of the AR model for arfit.m 
            'disc',100,...        %number of intervals of the Fourier domain 
            'norm',0,...          %normilize each trial
            'anom',0,...          %remove periodicity anom
            'test',0,...          %number of IAAFT data for sign-test
            'graph','y',...       %plot results y/n
            'vname',0);           %parse the name of the varaibles as cell array  
opts=parseArgs(varargin,opts);


%constant and parameters
m=opts.order;
dt=opts.dt;
[n,R]=size(x{1});
fn=1/2/dt;             %Niquist frequency
dth=pi/(opts.disc-1);  %interval of discretizzation of angular frequency
yax=[0 5];
M=length(x);

%variables names
if ~iscell(opts.vname)
    for i=1:M
        vname{i}=num2str(i);
    end
else
    for i=1:M
        vname{i}=opts.vname{i};
    end
end
    
max_iter=opts.max_iter;
tol=opts.tol;
anom=opts.anom;

theta=0:dth:pi;% discretizzation of fourier space (in angular frequency)
T=length(theta);
S=zeros(M,M,T);
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
f=zeros(M,M,T);

for j=1:M
FT{j}=zeros(T,n);
end

A=nan(M,m*M,R);
C=nan(M,M,R); 

  
for r=1:R
r/R
pause(.1)

    if opts.norm>0
        for j=1:M
           v(:,j)=red2(x{j}(:,r),'anom',anom);
        end

    else
        for j=1:M
           v(:,j)=x{j}(:,r);
        end
    end


%    for i=2:T
% 
%    s0=1/fr(i)/fourier_factor;
%      for j=1:M
%   [ft,period(i),scale,coi] = wavelet(v(:,j),dt,PAD,dj,s0,0,mother,param);
% 
%   FT{j}(i,:)=ft;
%      end
%    end
%    
%     for k=1:M
%         for j=1:M
%             f1(1,1,:)=mean(FT{k}.'.*conj(FT{j}).');
%             f(k,j,:)=f(k,j,:) + f1;
%         end
%     end
  
    WT = wave(v,'scale',sc(2:end),'dj',dj,'param',param,'dt',dt);
    f(:,:,2:end)=f(:,:,2:end)+WT;
    if strcmp(opts.extrap_method,'param')
    %fit parametric model for extrap
    [~, A(:,:,r), C(:,:,r)]=arfit(v,m,m);
    end
  
end

%extrapolate spectra at f = 0
if strcmp(opts.extrap_method,'param')
%using parametric method
  A(abs(A)==inf)=nan;
  S = AR_spectrum(nanmean(A,3),nanmean(C,3),fr*dt);S=S*dt;
  use=find(fr>=opts.fextrap & fr>fr(1)); 

  S(:,:,use)=f(:,:,use)/R/df;
else

use=find(fr>=opts.fextrap & fr>fr(1)); 
S(:,:,use)=f(:,:,use)/R/df;
S = spectral_matrix_extrap(f/R/df,fr,opts.fextrap);
end

% s11(:,1)=S(1,1,:);
% plot(fr,s11,'.-')
% pause
 [SIGMA4,H]=wilson(S,max_iter,tol);
 b=0;ylim=[0 0];
    for i=1:M
        for j=1:M
            b=b+1;
          if i==j  
          sii(:,1)=S(i,i,:);
           if opts.graph == 'y'
            subplot(M,M,b)
             plot(fr,sii);
            set(gca,'xlim',[0 fn])
            pause(.1)
           end
          else
              
%               if i<j
                  [SIGMA2,H2]=wilson(S([i j],[i j],:),max_iter,tol);
                  [F{i,j},Fxy{i,j},Fyx{i,j}] = Granger(SIGMA2,H2,S([i j],[i j],:)); 
%               else
%                   F{i,j}  = F{j,i};
%                   Fxy{i,j}= Fyx{j,i};
%                   Fyx{i,j}= Fxy{j,i};
%               end
          if opts.graph == 'y'
          subplot(M,M,b)
          plot(fr,Fyx{i,j},'k--');hold all
          plot(fr,F{i,j}+Fxy{i,j}+Fyx{i,j},'color',[0.4 0.4 0.4])
          set(gca,'xlim',[0 fn])
          ylim(1)=min(min(Fyx{i,j}),ylim(1));
          ylim(2)=max(max(Fyx{i,j}),ylim(2));
          end
          
          if M>2
            sub=find(1:M~=i & 1:M~=j);
          
            sub1=[i sub];
            sub2=[i j sub];
            [SIGMA3,G]=wilson(S(sub1,sub1,:),max_iter,tol);
            Fyx_z{i,j} = Granger_multicond(G,H(sub2,sub2,:),SIGMA3,SIGMA4(sub2,sub2));
          
            if opts.graph == 'y'
            plot(fr,Fyx_z{i,j},'k','linewidth',2)
            set(gca,'xlim',[0 fn])
            ylim(1)=min(min(Fyx_z{i,j}),ylim(1));
            ylim(2)=max(max(Fyx_z{i,j}),ylim(2));
            pause(.1)
            end
            
          end
            end

        end
    end
    
      if opts.graph=='y';b=0;
        for i=1:M
        for j=1:M
            b=b+1;
          if i~=j  

          subplot(M,M,b)
          set(gca,'ylim',ylim);
%   axis([0 3 -0.1 0.3])
          end
          
          if i==1
              subplot(M,M,b)
              title(['cause ' vname{j}])
          end
          
           if j==1
              subplot(M,M,b)
              ylabel(['effect ' vname{i}])
           end
          
        end
        end
      end
          

% if opts.test>0
%     
% h = waitbar(0,'IAAFT test, please wait ...');
% 
% for i=1:opts.test
%       waitbar(i/opts.test);
% %     surrogate=IAAFT(reshape(x,n*R,1));
% %     x_s=reshape(surrogate,n,R);
% %     Y=repmat(IAAFT(y(:,1)),1,R);
% %     [txy(:,i),tyx(:,i),txy_z(:,i),tyx_z(:,i)] =cond_spec_GC(x_s,y,z,'graph','n','disc',opts.disc);
% [txy(:,i),tyx(:,i),txy_z(:,i),tyx_z(:,i)] =cond_spec_GC(IAAFT(x),y,z,'graph','n','disc',opts.disc);
%     
% end
% close(h) 
% % varargout{1}=G12;
% % varargout{2}=G21;
% 
% 
%  if opts.test==1
%   subplot(subpl(2))
%   hold all
%   plot(fr,txy,'.-k')
%   plot(fr,tyx,'--k')
%   
%   subplot(subpl(3))
%   hold all
%   plot(fr,txy_z,'.-k')
%   plot(fr,tyx_z,'--k')
%  else
%   subplot(subpl(2))
%   hold all
%   plot(fr,prctile(txy',95),'.-k')
%   plot(fr,prctile(tyx',95),'--k')
%   
%   subplot(subpl(3))
%   hold all
%   plot(fr,prctile(txy_z',95),'.-k')
%   plot(fr,prctile(tyx_z',95),'--k')
%  end
% end

% 
function  S = spectral_matrix_extrap(S,fr,fextrap)
    if fextrap==0;fextrap=fr(2);end
gg=find(fr<=fextrap); %extrapolate after fextrap
hh=find(fr<fextrap/2 & fr>fextrap);%retain for interpolation
if isempty(hh)==1;hh=2;end
a=0;method='cubic';

for i=1:M
  for j=1:M
      
      if i==j
      yy(1,:)=S(i,j,hh);
      y1=fliplr(yy);
      y2=yy;
      x1=-fliplr(fr(hh));
      x2=fr(hh);
      xi=fr(gg);
      
    S(i,j,gg)=interp1([x1 x2],[y1 y2],xi,method);
    
      else
      yy1(1,:)=S(j,i,hh);
      yy2(1,:)=S(i,j,hh);
      y1=fliplr(yy1);
      y2=yy2;
      x1=-fliplr(fr(hh));
      x2=fr(hh);
      xi=fr(gg);
      
    S(i,j,gg)=interp1([x1 x2],real([y1 y2]),xi,method)+...
        sqrt(-1)*interp1([x1 0 x2],imag([y1 0 y2]),xi,method);
    
      end
%        a=a+1;
%   subplot(3,3,a)
%   ff(:,1)=S(i,j,:);semilogx(fr,[real(ff) imag(ff)],'.-')     
  end
  
end

end
end
