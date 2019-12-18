% compute interdependence and directional causality F, F12 and F21 (bivariate
% case only)
% from matrices:
% SIGMA : error prediction matrx (2,2)
% H     : transfer function matrx (2,2,:)
% S     : spectral matrx (2,2,:)

function [F,F12,F21] = Granger(SIGMA,H,S)

      
s11(1,:)=S(1,1,:);
s12(1,:)=S(1,2,:);
s21(1,:)=S(2,1,:);
s22(1,:)=S(2,2,:);

S11=SIGMA(1,1);
S12=SIGMA(1,2);
S21=SIGMA(2,1);
S22=SIGMA(2,2);


% % transform matrix to make the errors independent
% %you
% P1 = MVAR_norm_matrix(SIGMA);
% P2 = MVAR_norm_matrix(SIGMA([2 1],[2 1]));P2=P2';
% 
% 
% for i=1:length(H)
% H1(:,:,i)=H(:,:,i)*inv(P1);
% H2(:,:,i)=H(:,:,i)*inv(P2);
% end
% 
% H11(1,:)=H1(1,1,:);
% H22(1,:)=H2(2,2,:);
% 
% F21=log(s11./real(H11.*S11.*conj(H11)));
% F12=log(s22./real(H22.*S22.*conj(H22)));
% 
% F=log(real(H11.*S11.*conj(H11).*H22.*S22.*conj(H22))./real(s11.*s22-s12.*s21));

H11(1,:)=H(1,1,:);
H12(1,:)=H(1,2,:);
H21(1,:)=H(2,1,:);
H22(1,:)=H(2,2,:);

% F21=real(log(s11./(s11-(S22-S12^2/S11)*abs(H12).^2)));
% F12=real(log(s22./(s22-(S11-S21^2/S22)*abs(H21).^2)));

F21=log(abs(s11)./abs(s11-(S22-S12^2/S11)*abs(H12).^2));
F12=log(abs(s22)./abs(s22-(S11-S21^2/S22)*abs(H21).^2));

H11_t=H11+S12/S11*H12;
H22_t=H22+S21/S22*H21;

F=real(log(H11_t.*S11.*conj(H11_t).*H22_t.*S22.*conj(H22_t)./(s11.*s22-s12.*s21)));


