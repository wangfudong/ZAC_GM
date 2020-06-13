function p = normalize_unequal(q, group1, group2)
nn = size(q, 1);
n1 = size(group1, 2);
n2 = size(group2, 2);

X = zeros(n1, n2);
for ind = 1 : nn
    gid1 = group1(ind,:);
    gid2 = group2(ind,:);
    X(gid1, gid2) = q(ind);
end
X = bistocNormalize_slack(X, 1e-3);    

p = zeros(nn, 1);
for ind = 1 : nn
    gid1 = group1(ind,:);
    gid2 = group2(ind,:);
    p(ind) = X(gid1, gid2);
end

end