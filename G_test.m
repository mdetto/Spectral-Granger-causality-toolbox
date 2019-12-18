function [F1, F2, p_value]=G_test(x,y,m,p)
addpath ARFIT
%time domain G_test Bivariate Granger Causality Test with a priory
%defined order model m. Use F-stat to find p-value

%x and y two vector variables
%m order of the AR model
%p p-value for the F-test

N = length(y);
if length(x)~=N, error('Improper input args.'), end
if nargin<4, p=0.95; end

[~, ~, Cxy]=arfit([y,x], m, m);
[~, ~, Cy]=arfit(y, m, m);
[~, ~, Cx]=arfit(x, m, m);
 
RSS1=Cxy(1,1);
RSS0=Cy;
F1 = (N-2*m-1)*(RSS0-RSS1)/m/RSS0; %value of F-statistic
Gx_y=log(RSS0/RSS1);

RSS1=Cxy(2,2);
RSS0=Cx;

F2 = (N-2*m-1)*(RSS0-RSS1)/m/RSS0; %value of F-statistic
Gy_x=log(RSS0/RSS1);

 % look-up table for the inverse of F-statistic at specified p-value
sign = finv(1-p,m,N-2*m-1)
p_value= 1-fcdf(F2,m,N-2*m-1)
p_value= 1-fcdf(F1,m,N-2*m-1)


   if F2>F1
      disp('The second argument causes the first one.')
   else
      disp('The first argument causes the second one.')
   end

