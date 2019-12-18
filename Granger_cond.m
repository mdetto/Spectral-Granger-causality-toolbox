%from Chen (2006)
%Frequency decomposition of conditional Granger causality and application
%to multivariate neural field potential data

function F21_3=Granger_cond(S,H,SIGMA)


M=length(S);
for y=1:M

    
    H2=inv(H([1 3],[1 3],y))*[H(1,2,y);H(3,2,y)];
    
    SIGMA2=SIGMA([1 3],[1 3])+H2*[SIGMA(1,2) SIGMA(3,2)]+...
                 [SIGMA(1,2); SIGMA(3,2)]*H2'+...
                 SIGMA(2,2)*H2*H2';
             
    P=[1 0;-SIGMA2(1,2)/SIGMA2(1,1) 1];
    
    G(:,:,y)=H([1 3],[1 3],y)*inv(P);
    Sxx2(y)=SIGMA2(1,1);
    
    G2=[G(1,1,y) 0  G(1,2,y)
           0     1    0
    G(2,1,y)     0  G(2,2,y)];
    Q(:,:,y)=inv(G2)*H(:,:,y);
    
end

Qxx(1,:)=Q(1,1,:);
Sxx=SIGMA(1,1);

% F21_3=real(log(Sxx2./(Qxx.*Sxx.*conj(Qxx))));

F21_3=log(abs(Sxx2)./abs(Qxx.*Sxx.*conj(Qxx)));

