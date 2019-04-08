function w = weights(x)
% w = weights(x)
%
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S36144502417716

assert(size(x,2)==1);
n = numel(x);

w = zeros(size(x));
for j=1:n
    w(j) = 1/prod(x(j)-x(1:n~=j));
end