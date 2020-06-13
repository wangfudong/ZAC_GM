function X = IPFP_U(K, group1, group2)
x0 = ones(size(K,1), 1);
x0 = normalize(x0, group1);
X = IPFP(K, x0, group1, group2);
end