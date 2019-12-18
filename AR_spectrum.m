%function AR_spectrum form autoregressive coefficient
%
%  evaluete analitically the spectrum from the coefficient of an autoregressive model
%     v(k,:)' =A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)', 
%
%  where A=[A1 ... Ap] is the coefficient matrix, 
%  The vectors eta(k,:) are independent Gaussian noise
%  vectors with mean zero and covariance matrix C.
%  f is the adimentional frequency vector you want to evaluete the spectrum (for a unitary time step dt=1)

function [S,H,C] = AR_spectrum(A,C,f)

p=length(C); %number of variables 
m=length(A)/p;  %model order

for h=1:length(f)
    g=eye(p,p);
    for j=1:m
     g=g-A(:,(j-1)*p+1:j*p)*exp(-1i*2*pi*f(h)*j);
    end
    Hh=inv(g);
    S(:,:,h) = 1/pi * Hh*C*Hh';
    H(:,:,h) = pi^-0.5 * Hh;
%     S(:,:,h) =  Hh*C*Hh';
%     H(:,:,h) =  Hh;
end

%element on the diagonal are real
for j=1:p
S(j,j,:)=real(S(j,j,:));
end

% 
%  S = real(S);     
