%spectral matrix factorizzation S = psi*psi' , using Wilson algorytm (1972)
%
%
% S : nxn spectral matrix defined between [0 theta], theta = angular frequency
% e.g. S = [s11(theta) s12(theta)
%           s21(theta) s22(theta)]
%   s11 and s12 are the spectrum and co-spectrum between x1 and x2
%  
%
%      SIGMA : covariance matrix
%      H : spectral transfer function matrix (same dimention of S);
%      SIGMA=pi*A0*A0.';
%      H=psi*inv(real(A0));
%
%
%NOTE: the symbol " ' " denote matrix adjoint, 
%                 " .' " denote matrix transpose !!!
%
% Matteo Detto, PhD
% Biometeorology lab, Ecosystem Science Division
% Dept. of Environmental Science, Policy, and Management (ESPM)
% University of California, Berkeley, CA 94720
% Phone: 1-510-642 9048
% Fax: 1-510-643-5098 
% last update: 09 Oct 2009

function [SIGMA,H,B1,B2,B3]=wilson(S,varargin)

if nargin > 1
    max_iter=varargin{1};
    tol=varargin{2};
else
max_iter=50;
tol=1e-2;
end

SS=size(S);
m=SS(1);   %order od the system 
n=SS(3);   %length of the frequency domain
dth=pi/(n-1);
theta=-2*pi:dth:2*pi;
t=[0:dth:pi];
N=length(theta);
M=length(t);


% build hermitian matrix (f(theta)=f(-theta).') 
% and periodic between [-2pi 2pi] with period 2pi
b=0;
for j=ceil(N/2):ceil(N*3/4)
    b=b+1;
    
    for k=1:m
        for l=1:m  
           f(l,k,j)=S(k,l,b);
        end
    end

    j2=j-floor(N/2);
    for k=1:m
        for l=1:m  
          f(l,k,j2)=S(k,l,b);
        end
    end

end
%------------
for j=ceil(N/4):ceil(N/2)
    
    for k=1:m
        for l=1:m  
          f(l,k,j)=S(l,k,b);
        end
    end

    j2=j+floor(N/2);
    for k=1:m
        for l=1:m  
          f(l,k,j2)=S(l,k,b);
        end
    end
   b=b-1;
end

% ceil(N/2)
% se=find(theta==0)
% f(:,:,se)
% pause
% f(:,:,ceil(N/2))=real(f(:,:,ceil(N/2)));

% f(:,:,ceil(N/2)+100)
% f(:,:,1+100)
% pause
% %check if the matrix is hermitian and periodic
% for qq=1:100
% GG(:,:)=f(:,:,ceil(N/2)+qq);
% FF(:,:)=f(:,:,ceil(N/2)-qq).';
% fdf=abs(GG(1,2)-FF(1,2));
% if real(fdf)>0 | imag(fdf)>0
%     'matrix not-hemitian'
%     pause
% end
% 
% end
% 
% %plot co-spectral matrix
% figure(2)
% clf
% subplot(221)
% AWE(:,1)=f(1,1,:);
% plot(theta,AWE)
% subplot(224)
% AWE(:,1)=f(2,2,:);
% plot(theta,AWE)
% subplot(222)
% AWE(:,1)=real(f(1,2,:));AWE(:,2)=imag(f(1,2,:));
% plot(theta,AWE)
% subplot(223)
% AWE(:,1)=real(f(2,1,:));AWE(:,2)=imag(f(2,1,:));
% plot(theta,AWE)
% pause

use=find(theta>=-pi & theta<=pi);

%check if the integral of the log of det(f)>-inf
for w=1:length(use)
    Y(w)=det(f(:,:,use(w)));
end

if real(trapz(theta(use),log(Y)))==-inf
    'warning: the integral of the log of det(f)=-inf'
end
%-----------------------

% %first guess
g0=real(1/2/pi*trapz(theta(use),f(:,:,use),3));
A0 =inv(chol(inv(g0)));
psi1=repmat(A0,[1 1 N]);
      
use=find(theta>=-pi & theta<=pi);
%iteration loop
iter=0;
%G0_old=eye(m,m);
%G0=2*eye(m,m);
R=zeros(m,m);


AWE(:,1)=f(1,1,randperm(N));
AWE(:,2)=f(1,1,:);
rsq=corrcoef(AWE(:,1),AWE(:,2));

%while (iter<=max_iter & abs(1-abs(det(G0_old)/det(G0)))>tol)
 while (iter<=max_iter && 1-rsq(2,1)^2>tol)
    iter=iter+1;
    psi0=psi1; 


% 1-abs(det(G0_old)/det(G0))
rsq=real(corrcoef(AWE(:,1),AWE(:,2)));
% 1-rsq(2,1).^2
% pause

%G0_old=G0;
G0=real(A0*A0.');
g=zeros(m,m,N);

 if iter==max_iter
   display( ['max iteration reached:   ' num2str(max_iter)] )
   break
 end

   %  g-function with periodic conditions     
     for j=ceil(N/4):ceil(N/2)
%      g(:,:,j)=inv(psi0(:,:,j))*f(:,:,j)*inv(psi0(:,:,j))'+eye(m);
       g(:,:,j)=psi0(:,:,j)\f(:,:,j)/psi0(:,:,j)'+eye(m);
       g(:,:,N-j+1)=g(:,:,j).';
     end

       g(:,:,ceil(N*3/4)+1:N)=g(:,:,ceil(N/4)+1:ceil(N/2));
       g(:,:,1:ceil(N/4)-1)  =g(:,:,ceil(N/2):ceil(N*3/4)-1);

    %theta loop, from 0 to pi
    b=0;
for j=ceil(N/2):ceil(N*3/4)
   %integrand for h-function calculation
   integrand(:,:,1)=2*(g(:,:,j-1)-g(:,:,j+1))/dth; 
   integrand(:,:,2:M)=0;  
   for k=2:M
      integrand(:,:,k)=(g(:,:,j-k+1)-g(:,:,j+k-1))*cot(t(k)/2);
   end
 
 % h-function, solve integral with trapezodial approx. 
   h=1/2/pi*trapz(t,integrand,3);

  % g_plus-function
  g_plus(:,:)=1/2*(g(:,:,j)+1i*h(:,:));
  
  %create a matrix R so psi1 upper triangular, which satisfy R + R' = 0
  if j==ceil(N/2)
      for k=1:m
          for l=1:m
              if k<l
   R(k,l)=conj(g_plus(l,k));
   R(l,k)=-g_plus(l,k);
              end
          end
      end
  end

  psi1(:,:,j)=psi0(:,:,j)*(g_plus + R);
  psi1(:,:,ceil(N/2)-b)=conj(psi1(:,:,j));
  b=b+1;
end

       %make psi1 periodic
       psi1(:,:,1:ceil(N/4)-1)  =psi1(:,:,ceil(N/2):ceil(N*3/4)-1);
       psi1(:,:,ceil(N*3/4)+1:N)=psi1(:,:,ceil(N/4)+1:ceil(N/2));

% %uncomment this block if you want to see how the solution approach to the real
for y=1:N
FDF(:,:,y)=psi1(:,:,y)*psi1(:,:,y)';
end
%        
% figure(1)
% clf
% subplot(221)
AWE(:,1)=FDF(1,1,:);
AWE(:,2)=f(1,1,:);
% plot(theta(use),AWE(use,:))
% legend('model','obs')
% pause(0.1)
% subplot(224)
% AWE(:,1)=FDF(2,2,:);
% AWE(:,2)=f(2,2,:);
% plot(theta,AWE)
% 
% subplot(222)
% AWE(:,1)=real(FDF(1,2,:));AWE(:,2)=real(f(1,2,:));
% plot(theta,AWE)
% hold on
% AWE(:,1)=imag(FDF(1,2,:));AWE(:,2)=imag(f(1,2,:));
% plot(theta,AWE)
% subplot(223)
% AWE(:,1)=real(FDF(2,1,:));AWE(:,2)=real(f(2,1,:));
% plot(theta,AWE)
% hold on
% AWE(:,1)=imag(FDF(2,1,:));AWE(:,2)=imag(f(2,1,:));
% plot(theta,AWE)
% pause

A0=1/2/pi*trapz(theta(use),psi1(:,:,use),3);
%A0(:,:,1)=psi1(:,:,ceil(N/2));

end


for k=1:m
    for l=1:m
     PSI(1,:)=psi1(k,l,use);
     A1(k,l)=1/2/pi*trapz(theta(use),PSI.*exp(-1i*1*theta(use)));
     A2(k,l)=1/2/pi*trapz(theta(use),PSI.*exp(-1i*2*theta(use)));
     A3(k,l)=1/2/pi*trapz(theta(use),PSI.*exp(-1i*3*theta(use)));
    
    end
end
B1=real(A1/A0);
B2=real(A2/A0);
B3=real(A3/A0);


%--------------------------
SIGMA=pi*real(A0*A0.');

for u=ceil(N/2):ceil(N*3/4)
 H(:,:,u-ceil(N/2)+1)=pi^-0.5*psi1(:,:,u)/real(A0);
end

