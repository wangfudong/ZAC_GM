function X = PSM(K, group1, group2, iter)
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
%   X: steady state distribution of PSM
%
% History
%   create  -  A. Egozi
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 01-03-2015
    
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
    
    for i = 1 : iter
        p_ = p;
        q = AA * p;
        
        if rem(i, 2) == 1
            p = normalize(q, group1);
        else
            p = normalize(q, group2);
        end
        p_ = p_ + eps;
        
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
            fprintf('\nPSM iterator = %d', i);
        end

    end
    X = p;         
end
