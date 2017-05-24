function [ out ] = cnormdpp( sz, k, bordergap )
%CNORMDPP   k-sparse DPP pattern using eigenvalues of iid c-normal matrix.
    if nargin < 3
        bordergap = 1.1;
    end

    if size(sz) ~= 2
        error('CNORMDPP only creates 2D activation maps');
    end

   M = randn(k) + 1i*randn(k);
   s = eig(M);
   
   x = real(s);
   x = x/(2*norm(x, 'inf')*bordergap) + 0.5;

   y = imag(s);
   y = y/(2*norm(y, 'inf')*bordergap) + 0.5;

   x = ceil(x * sz(1));
   y = ceil(y * sz(2));
   out = zeros(sz);
   out(sub2ind(sz, x, y)) = 1;
end

