function nF = nF_by_kmeans(X,Y,SC_or_dis,rota)

if SC_or_dis == 0
    M = M_shape(X,Y,1/8,2,rota);
else
    M = M_points(X,Y);
end
M = M/max(M(:));

M_min = zeros(length(X(:,1)),1);
for i = 1:length(X(:,1))
    M_min(i) = min(M(i,:));
end
% [val,ind1] = vl_kmeans(M_min',2);
[ind1,val] = kmeans(M_min,2);
val1 = find(val==min(val));

M_min = zeros(length(Y(:,1)),1);
for i = 1:length(Y(:,1))
    M_min(i) = min(M(:,i));
end

% M0 = asgHun(-M);
% M_min = M(M0==1);

% [val,ind2] = vl_kmeans(M_min',2);
[ind2,val] = kmeans(M_min,2);
val2 = find(val==min(val));

nF=min(sum(ind1==val1),sum(ind2==val2));

