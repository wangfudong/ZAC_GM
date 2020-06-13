function X = PSM_EX(K, group1, group2, iter)
% Probabilistic to Spectral Matching.
%
% Reference
%   A. Egozi, Y. Keller, and H. Guterman. "A Probabilistic Approach to
%   Spectral Graph Matching", IEEE Trans. Pattern Anal. Mach. Intell.,
%   35(1):18-27, 2013.
%
% Math
%   This algorithm is to obtain the optimal X for the following problem
%     max_x   x' * K * x
%     s.t.    A * x <= 1
%
% Input
%   K:      [m * m] affinity matrix of m candicate matches
%   group1: conflicting match groups in domain 1 (size(K,1) x nGroup1)
%   group2: conflicting match groups in domain 2 (size(K,1) x nGroup2)
%                 
%       e.g. find(group1(:,3)) represents the third goup of matches  
%                               sharing the same feature in domain1   
%
% Output
%   X: steady state distribution of PSM_EX
%
% History
%   create  -  Tao WANG (twang@bjtu.edu.cn), 05-11-2015



    mm = size(K, 1);
    ro = 1e-6;
    if nargin < 4 || iter <= 0
        iter = 1000;
    end
    
    p = ones(mm, 1);
    if issparse(K)
        p = sparse(p);
    end
%     p = normalize(p, group1);
    AA = K;
    
%     % bistoc normalize affinity matrix by edge
%     AA = bistocNormalizeAffinity(AA, group1, group2);
   
    for i = 1 : iter
        
        p_ = p;
        q = AA * p;
        
        % impose 1-1 match constraints
        p = bistocNormalizeSolution(q, group1, group2);
%         if rem(i, 2) == 1
%             p = normalize(q, group1);
%         else
%             p = normalize(q, group2);
%         end
        p_ = p_ + eps;
        
        
        % refine affinity matrix
        if issparse(K)
            AA = AA';
            for col = 1:size(AA,1)
                AA(:,col) = AA(:,col) *(p(col) / p_(col));
            end
            AA = AA';
        else
            AA = AA .* ((p./ p_) * (ones(1, mm)));
        end

        if sqrt(sum((p - p_).^2))/(mm) < ro
            break;   
        end 
        
        if mod(i,100) == 0
            fprintf('\nPSM_EX iterator = %d', i);
        end

    end
    X = p;     
    
end

function K = bistocNormalizeAffinity(K, group1, group2)
%% generate nodes and matches info
[ind_match, ind_grp1] = find(group1);
grp1 = [ind_match, ind_grp1];
grp1 = sortrows(grp1, 1);
grp1 = grp1(:,2);
N1 = max(grp1);

[ind_match, ind_grp2] = find(group2);
grp2 = [ind_match, ind_grp2];
grp2 = sortrows(grp2, 1);
grp2 = grp2(:,2);
N2 = max(grp2);

match = zeros(N1, N2);
for k = 1:size(grp1, 1)
    match(grp1(k), grp2(k)) = k;
end

%% convert affinity matrix to edge similarity matrix
edges1 = zeros(N1 * N1, 4);     %[index, tail_node, head_node, weight]
edges1(:,1) = 1:N1 * N1;
edges1(:,2) = ceil(edges1(:,1) / N1);
edges1(:,3) = mod(edges1(:,1) - 1, N1) + 1;
    
edges2 = zeros(N2 * N2, 4);     %[index, tail_node, head_node, weight]
edges2(:,1) = 1:N2 * N2;
edges2(:,2) = ceil(edges2(:,1) / N2);
edges2(:,3) = mod(edges2(:,1) - 1, N2) + 1;
    
S = zeros(N1 * N1, N2 * N2);
[mi, mj, v] = find(K);
for k = 1:size(mi,1)
    ni = grp1(mi(k));   % node ni in G1 match to node ni_ in G2
    ni_ = grp2(mi(k));
    nj = grp1(mj(k));   % node nj in G1 match to node nj_ in G2
    nj_ = grp2(mj(k));
        
    S((ni - 1) * N1 + nj, (ni_ - 1) * N2 + nj_) = v(k);
end
edges1(:,4) = sum(S, 2);
edges2(:,4) = sum(S, 1);
    
%% remove empty edges 
S = S(:, any(S, 1));
S = S(any(S, 2), :);
    
edges1 = edges1(any(edges1(:,4), 2), :);
edges2 = edges2(any(edges2(:,4), 2), :);
    
%% bistoc normalize edge similarity
tolC = 1e-3;
S = bistocNormalize_slack(S, tolC);
    
%% convert similarity matrix back to affinity matrix
K = zeros(size(K));
[e1, e2, v] = find(S);
for k = 1:size(e1,1)
    ni = edges1(e1(k), 2);
    nj = edges1(e1(k), 3);
    ni_ = edges2(e2(k), 2);
    nj_ = edges2(e2(k), 3);
    ind_i = match(ni, ni_);
    ind_j = match(nj, nj_);
    if (ind_i <= 0 || ind_j <= 0)
        printf('\r\n ni = %d, nj = %d, ni_ = %d, nj_ = %d', ni, nj, ni_, nj_);
        error(['ind_i = ' num2str(ind_i) ', ind_j = ' num2str(ind_j) ', error!']);
    end
    K(ind_i, ind_j) = v(k);
end 

end

function p = bistocNormalizeSolution(p, group1, group2)
%% generate nodes and matches info
[ind_match, ind_grp1] = find(group1);
grp1 = [ind_match, ind_grp1];
grp1 = sortrows(grp1, 1);
grp1 = grp1(:,2);
N1 = max(ind_grp1);

[ind_match, ind_grp2] = find(group2);
grp2 = [ind_match, ind_grp2];
grp2 = sortrows(grp2, 1);
grp2 = grp2(:,2);
N2 = max(ind_grp2);

%% convert identity vector to assignment matrix
X = zeros(N1, N2);
for i = 1:size(p)
    X(grp1(i), grp2(i)) = p(i);
end

%% bistoc normalize assignment matrix
tolC = 1e-3;
X = bistocNormalize_slack(X, tolC);

%% convert assignment matrix back to identity vector
p = zeros(size(group1,1), 1);
for i = 1:size(p)
    p(i) = X(grp1(i), grp2(i));
end

end
