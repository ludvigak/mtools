function fi = direct(x, w, f, xi)
% fi = bclag_interp_direct(x, w, f, xi)
%
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715


n = numel(x);
N = numel(xi);
dim = size(f, 2);

exact = zeros(N, 1);
numer = zeros(N, dim);
denom = zeros(N, dim);

for j=1:n
    xdiff = xi-x(j);
    temp = w(j)./xdiff;
    numer = numer + temp*f(j,:);
    denom = denom + temp;
    exact(xdiff==0) = j;
end
fi = numer./denom;
jj = find(exact); 
fi(jj,:) = f(exact(jj),:);