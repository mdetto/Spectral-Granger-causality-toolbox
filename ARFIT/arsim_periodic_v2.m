function [v]=arsim_periodic(w,omega,Ampl,Fi,C,n)
%ARSIM	Simulation of AR process with one periodic
%
%  v=ARSIM(w,A,C,n) simulates n time steps of the AR(p) process
%
%     v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)', 
%
%  Modified 13-March-2010
%  Author: Matteo Detto

  m=length(w);
  w       = w(:)';                      % force w to be row vector

  % Compute Cholesky factor of covariance matrix C
  [R, err]= chol(C);                    % R is upper triangular
  if err ~= 0
    error('Covariance matrix not positive definite.')
  end
    
  %define error predictions
  for i=1:m
  a=linspace(Ampl(i,1),Ampl(i,2),n);
  f=linspace(Fi(i,1),Fi(i,2),n);
  randamp = a(randperm(n))';
  randFi  = f(randperm(n))';
  prdc=randamp.*sin(omega*[1:n]' + randFi);

  Wien(:,i)=prdc;
  end
  
  % Get ndisc+n independent Gaussian pseudo-random vectors with 
  % covariance matrix C=R'*R
  randvec = Wien*R;

  % Add intercept vector and periodicity to random vectors
  v = randvec + ones(n,1)*w ;

