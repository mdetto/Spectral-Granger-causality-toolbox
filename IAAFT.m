% generate Foutier Trasform surrogate with iterative method
% iterative amplitude adjusted Fourier transform

function [iaaft,index] = IAAFT(template,varargin)

%default options
opts=struct('mother','Morlet',...  
            's0',-1, ...    
            'tol',1e-3, ...    
            'param',15, ...    
            'max_iter',100, ...    
            'dj',0.25, ... 
            'dt',1,...
            'order',2,...
            'disc',2);
opts=parseArgs(varargin,opts);

tol=opts.tol;
maxiter=opts.max_iter;

[n,R]=size(template);


%param&const
dj=0.15;
k0=6;
Ff=4*pi/(k0+sqrt(2+k0^2));
s0=2/Ff;
coi=((n+1)/2-1)/sqrt(2); 
J1=floor((log(coi/s0)/log(2))/dj);
scale(:,1) = s0*2.^((0:J1)*dj);
npuls  = floor((n-1)/2);
K(:,1)  = 2*pi/n*[ 0:npuls  (npuls-n+1):-1 ];
H=exp(-(K*scale'-k0).^2);

size(H)
iaaft=zeros(n,R); 
for k=1:R

FT = fft(template(:,k));
sortedValues=sort(template(:,k));
fourierCoeff=abs(FT);

surrogate=template(randperm(n),k);

%inizialize variables
FS=FT.*conj(FT)/n;
E0(:,1)=(FS.'*H)./sum(H);
 
iter=0;p=0;
r=0;

while 1-r^2>tol && iter<maxiter && r<1

spectrum = ifft(surrogate);              % surrogate is the surrogate time series.
phase    = angle(spectrum);              % Angle is the Matlab function to calculate the phase from a complex number.
spectrum = fourierCoeff .* exp(1i*phase); % fourierCoeff are the magnitudes of the Fourier coefficients of the template.
surrogate = fft(spectrum);

[~, index]   = sort(real(surrogate));  % We need only the indices. The first value of index points to the highest value, etc.
surrogate(index) = sortedValues;     % sortedValues is the vector with the sorted values from the template.

FT1=fft(surrogate);
FS1=FT1.*conj(FT1)/n;
E1(:,1)=(FS1.'*H)./sum(H);

iter=iter+1;
r=corr(log(E0),log(E1));

% figure(1);clf
% plot(log(E0));hold all
% plot(log(E1))
% 1-r
% pause

if iter==maxiter
   
    surrogate=template(randperm(n),k);
    iter=1;
    p=p+1;
    disp(['max iterations reached:  ' num2str(maxiter) ' : ' num2str(p)])
end

end

iaaft(:,k)=surrogate;
end


    
