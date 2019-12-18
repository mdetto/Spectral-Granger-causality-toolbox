% compute the multi-variate conditional spectral G-causality
% extending Ding et al., 2005
% G: bivariate spectral transfer function with error matrix SIGMA3
% H: multivariate spectral transfer function with error matrix SIGMA4


function Fyx_z = Granger_multicond(G,H,SIGMA3,SIGMA4)

sz=size(H);
M=sz(1);
T=sz(3);

%transform matrix to make the errors independent
P3 = MVAR_norm_matrix(SIGMA3);
P4 = MVAR_norm_matrix(SIGMA4);

for i=1:T
G1(:,:,i)=G(:,:,i)*inv(P3);
H1(:,:,i)=H(:,:,i)*inv(P4);
end


 %expand G matrix
 sub=find(1:M~=2);
 G2=zeros(sz);
 G2(2,2,:)=ones(1,1,T);
 G2(sub,sub,:)=G1;
 
 for y=1:T
         Q(:,:,y)=inv(G2(:,:,y))*H1(:,:,y);
 end
 
Qxx(1,:)=Q(1,1,:);
Sxx=SIGMA4(1,1);
S3=SIGMA3(1,1);

Fyx_z=log(S3./abs(Qxx.*Sxx.*conj(Qxx)));
 