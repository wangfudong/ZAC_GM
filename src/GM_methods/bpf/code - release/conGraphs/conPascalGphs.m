function [K, group1, group2, KP, KQ, gphs, X_GT] = conPascalGphs(fmat, outliers)
load(fmat);
features1(:,3) = features1(:,9);     % for convenience, use the 3 column to store the orientation of the contour
features2(:,3) = features2(:,9);

%% generate feature points
gts = size(gTruth, 2);
gphs{1}.Pts = features1(1:gts, 1:3);
gphs{2}.Pts = features2(1:gts, 1:3);

% randomly add outliers
if outliers > 0
    k1 = min(outliers, size(features1, 1) - gts);
    k2 = min(outliers, size(features2, 1) - gts);
    p1 = randperm(size(features1, 1) - gts, k1) + gts;
    p2 = randperm(size(features2, 1) - gts, k2) + gts;
    gphs{1}.Pts = [gphs{1}.Pts; features1(p1, 1:3)];
    gphs{2}.Pts = [gphs{2}.Pts; features2(p2, 1:3)];
end

% exchange x and y coordinate
gphs{1}.Pts(:, 1:2) = [gphs{1}.Pts(:,2), gphs{1}.Pts(:,1)];
gphs{2}.Pts(:, 1:2) = [gphs{2}.Pts(:,2), gphs{2}.Pts(:,1)];

% randomly reorder points
n1 = size(gphs{1}.Pts, 1);
n2 = size(gphs{2}.Pts, 1);
ord1 = randperm(n1);
ord2 = randperm(n2);
gphs{1}.Pts = gphs{1}.Pts(ord1, :);
gphs{2}.Pts = gphs{2}.Pts(ord2, :);

% ground truth assignment
X_GT = zeros(n1, n2);
X_GT(1:gts, 1:gts) = eye(gts);
X_GT = X_GT(ord1, ord2);


%% build graphs
A1 = DelaunayTriGraph(gphs{1}.Pts);
A2 = DelaunayTriGraph(gphs{2}.Pts);

% edges
[p1, q1, ~] = find(A1);
[p2, q2, ~] = find(A2);
gphs{1}.Eg = [p1'; q1'];
gphs{2}.Eg = [p2'; q2'];

% matrix G and H used in FGM algorithm
[p1,q1,~] = find(tril(A1));
[p2,q2,~] = find(tril(A2));
m1 = size(p1,1);    
m2 = size(p2,1);

gphs{1}.G = zeros(n1, m1);
gphs{2}.G = zeros(n2, m2);
for i = 1:m1
    gphs{1}.G(p1(i), i) = 1;    
    gphs{1}.G(q1(i), i) = 1;
end
for i = 1:m2
    gphs{2}.G(p2(i), i) = 1;
    gphs{2}.G(q2(i), i) = 1;
end
gphs{1}.H = [gphs{1}.G eye(n1)];
gphs{2}.H = [gphs{2}.G eye(n2)];


%% compute node affinity
%K = zeros(n1 * n2, n1 * n2);
K = sparse(n1 * n2);
KP = zeros(n1, n2);
for i = 1:n1
    for k = 1:n2
        ori_i = gphs{1}.Pts(i, 3);
        ori_k = gphs{2}.Pts(k, 3);
        KP(i,k) = exp(-abs(ori_i - ori_k));
    %    KP(i,k) = exp(-abs(ori_i - ori_k)/pi);
    
        match_ind = (k - 1) * n1 + i;
        K(match_ind, match_ind) = KP(i,k);
    end
end


%% compute edge affinity
S1 = nanstd(gphs{1}.Pts(:, 1:2));
S2 = nanstd(gphs{2}.Pts(:, 1:2));
KQ = zeros(m1, m2);
for i = 1:m1
    cor_1p = gphs{1}.Pts(p1(i), 1:2);
    cor_1q = gphs{1}.Pts(q1(i), 1:2);
%    y_dist = abs(cor_1p(1) - cor_1q(1));
%    x_dist = abs(cor_1p(2) - cor_1q(2));
    y_dist = abs(cor_1p(1) - cor_1q(1)) / S1(1);
    x_dist = abs(cor_1p(2) - cor_1q(2)) / S1(2);
    dst1 = sqrt(y_dist * y_dist + x_dist * x_dist);
    ang1 = atan(y_dist / (x_dist + eps));

    u1 = p1(i);
    v1 = q1(i);
    for k = 1:m2
        cor_2p = gphs{2}.Pts(p2(k), 1:2);
        cor_2q = gphs{2}.Pts(q2(k), 1:2);
    %    y_dist = abs(cor_2p(1) - cor_2q(1));
    %    x_dist = abs(cor_2p(2) - cor_2q(2));
        y_dist = abs(cor_2p(1) - cor_2q(1)) / S2(1);
        x_dist = abs(cor_2p(2) - cor_2q(2)) / S2(2);
        dst2 = sqrt(y_dist * y_dist + x_dist * x_dist);
        ang2 = atan(y_dist / (x_dist + eps));
        
        Dst = abs(dst1 - dst2);
        Ang = abs(abs(ang1) - abs(ang2));
    %    Dst = abs(dst1 - dst2) / (max(dst1, dst2) + eps);
    %    Ang = abs(ang1 - ang2) / pi;
        KQ(i,k) = exp(-(Dst + Ang) / 2);
        
        u2 = p2(k);
        v2 = q2(k);
        
        ind_u1u2 = (u2 - 1) * n1 + u1;      % u1 --> u2
        ind_u1v2 = (v2 - 1) * n1 + u1;      % u1 --> v2 
        ind_v1v2 = (v2 - 1) * n1 + v1;      % v1 --> v2
        ind_v1u2 = (u2 - 1) * n1 + v1;      % v1 --> u2
        K(ind_u1u2, ind_v1v2) = KQ(i,k);
        K(ind_v1v2, ind_u1u2) = KQ(i,k);
        K(ind_u1v2, ind_v1u2) = KQ(i,k);
        K(ind_v1u2, ind_u1v2) = KQ(i,k);
    end
end


%% group info
Ct = ones(size(KP));
[ind, ind_m] = find(Ct);
%group1 = zeros(size(ind, 1), n1);
%group2 = zeros(size(ind, 1), n2);
group1 = sparse(size(ind, 1), n1);
group2 = sparse(size(ind, 1), n2);
for k = 1:size(ind, 1)
    group1(k, ind(k)) = 1;
    group2(k, ind_m(k)) = 1;
end

end