function M = Dis_sift(DE1,DE2)
m = size(DE1,1);
n = size(DE2,1);
M = zeros(m,n);

for i = 1:128
    M = bsxfun(@minus, DE1(:,i), DE2(:,i)').^2 + M;
end
M = M.^(0.5);
M = M/max(M(:));