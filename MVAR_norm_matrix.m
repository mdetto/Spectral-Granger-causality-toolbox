% find the matrix P to Trasform a MVAR of order m with cross-correlated noise 
% to a MVAR with uncorrelated noise. The normalizzation is made left-multiplying
% by P MVAR in the form B(L)Y = E ;
% this function is needed for multi_spec_GC, but it can be used with
% spec_GC. The trick of the normalizzation permit to separate the pure directional interactions 

function P = MVAR_norm_matrix(C)

m=length(C);

P=eye(m);
 for i=1:m-1
   
     Pi=eye(m);
     
     for j=i+1:m
         Pi(j,i)=-C(j,i)/C(i,i);
     end
  
%    Pi
%    pause
     for p=i+1:m
         for q=i+1:m
             C(p,q)=C(p,q)-C(p,i)/C(i,i)*C(i,q);
         end
     end
     
     C(i,i+1:m)=0;
     C(i+1:m,i)=0;
     
%      C
%      Pi
%      pause
     P=P*Pi;
 end
 
     
     



