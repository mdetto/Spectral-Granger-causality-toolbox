%compute discrete Fourier transform (not FFT)
% v is a vector of data
% theta is the angular frequency
function FT = fourier(v,theta)
  
  vn(:,1)=v-mean(v); %remove mean
  N=length(v);
  T=length(theta);
  n(1,:)=0:N-1;
  
  FT=zeros(T,1);
  for j=1:T
      
      FT(j)=exp(-theta(j)*1i*n)*vn;
  end