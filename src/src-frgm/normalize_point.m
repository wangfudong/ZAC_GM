function [Xn,max_bound,ux] = normalize_point(X,s)
if nargin < 2 || isempty(s)
    s = 1;
end
[m,n] = size(X);
bound = zeros(1,n);
for i = 1:n
    bound(i) = max(max(X(:,i)))-min(min(X(:,i)));
end
max_bound = max(bound);
ux = mean(X);
X = X - repmat(mean(X),m,1);
Xn = s*X/max_bound;