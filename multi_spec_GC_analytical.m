% Compute spectral Granger causality for the multi-variate analitically 
% from the coefficient matrix A of an AR model and relate error matrix C
% parseArgs.m is used to pass additional options
%
% [Fyx,Fxy,Fyx_z,varargout] = multi_spec_GC_analytical(A,C,varargin)
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

function [Fyx,Fxy,Fyx_z,varargout] = multi_spec_GC_analytical(A,C,varargin)

%default options
opts=struct('dt',1,...            %temporal step
            'max_iter',20, ...    %max iterations for Wilson algorithim
            'tol',1e-6, ...       %tolerance for Wilson algorithim
            'disc',100,...        %number of intervals of the Fourier domain 
            'graph','y',...       %plot results y/n
            'subpl',321,...
            'vname',[{'x'},{'y'},{'z'}]);%parse the name of the varaibles as cell array 
opts=parseArgs(varargin,opts);


%constant and parameters
dt=opts.dt;
M=length(C);
fn=1/2/dt;             %Niquist frequency
dth=pi/(opts.disc-1);  %interval of discretizzation of angular frequency
yax=[0 5];

max_iter=opts.max_iter;
tol=opts.tol;

theta=0:dth:pi;% discretizzation of fourier space (in angular frequency)
T=length(theta);
S=zeros(2,2,T);
fr=theta/2/pi/dt;

[S,H,SIGMA4] = AR_spectrum(A,C,fr);
S=S*dt;
H=conj(H)*sqrt(dt);
% [SIGMA4,H2]=wilson(S,max_iter,tol);
% figure(1);clf
% b=0;
% for i=1:3
%     for j=1:3
%     b=b+1;
%     subplot(3,3,b)
% s11(:,1)=imag(H(i,j,:));
% s11(:,2)=imag(H(i,j,:));
% plot(fr,s11,'.-')
% pause(.1)
% end
% end
% pause
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
    
        b=0;
        for i=1:M
        for j=1:M
            b=b+1;
          if i~=j  

          subplot(M,M,b)
          set(gca,'ylim',ylim);
          end
          
          if i==1
              subplot(M,M,b)
              title(['cause ' opts.vname{j}])
          end
          
           if j==1
              subplot(M,M,b)
              ylabel(['effect ' opts.vname{i}])
           end
          
        end
        end
          