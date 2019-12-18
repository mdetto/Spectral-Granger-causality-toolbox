function f = daub(k)

% Length of R(z), which is also the order of the filter
r = 2*k-1;

% Let's compute the coefficients of (1 + z^-1)^k * (1 + z)^k
% Remember that polynomial multiplication is equivalent to convolution

p = [1 1]; 	%  This can represent both (1 + z^-1) and (1 + z)

for n=1:k-1
	p = conv(p, [1 1]);
end
pz = conv(p, p);

% Now we want to build a convolution matrix A so that AR , with 
%  R = [r(k) r(k-1) ... r(0) ... r(k)]', gives us the coefficients of
%  the even powers of P(z).

% Time-reverse (not necessary due to symmetry) and zero-pad pz
pz = [zeros(1,r) pz zeros(1,r)];
l = length(pz)

A = zeros(r,r);
for n=1:r
	A(n,:) = pz(l-2*n-r+1:l-2*n);
end

% Now we solve AR = B, where B is zero for all odd powers except for z^0

B = zeros(r,1); B((r+1)/2) = 1;
R = A\B;
