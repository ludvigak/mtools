function B = matrix_nochecks(x, w, xi)
% B = bclag_interp_matrix_nochecks(x, w, xi)
% Same as bclag_interp.matrix without checks to see if any points in x and xi coincide
%
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715

assert(size(x,2)==1)
assert(size(xi,2)==1)
assert(all(size(x)==size(w)));

n = numel(x);
N = numel(xi);

B = zeros(N,n);
denom = zeros(N,1);

for j=1:n
    xdiff = xi-x(j);
    temp = w(j)./xdiff;
    B(:,j) = temp;
    denom = denom + temp;
end

B = bsxfun(@rdivide,B,denom);
