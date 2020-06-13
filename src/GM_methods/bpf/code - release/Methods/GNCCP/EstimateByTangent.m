function x = EstimateByTangent(x_store, sLen_store, bRefine)
%% estimate next solution x according to previous K solutions
% x_store:  [N * K] matrix
% sLen_store: K vector of step length

if issparse(x_store)
    x_shift = sparse(size(x_store,1), size(x_store,2));
else
    x_shift = zeros(size(x_store));
end

%% estimate shift of next solution
k = size(x_store, 2);
weight = 1:k-1;

rho = 10 * k;
weight = exp(-(k-1-weight).^2/rho);

for i = 1:k-1
    x_shift(:,i) = (x_store(:,i+1) - x_store(:,i)) / sLen_store(i);
    x_shift(:,k) = x_shift(:,k) + weight(i) * x_shift(:,i);
end
x_shift(:,k) = x_shift(:,k) / sum(weight);

x = x_store(:,k) + x_shift(:,k) * sLen_store(k);

if bRefine
    x(x<0)=0;
    x(x>1)=1;
end

end