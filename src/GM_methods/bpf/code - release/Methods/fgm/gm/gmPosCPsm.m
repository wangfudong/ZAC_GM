function X = gmPosCPsm(K, X0, par)
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
% Remark
%   nn = n1 x n2
%
% Input
%   K       -  affinity matrix, nn x nn (sparse)
%   ns      -  #nodes, 1 x 2
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg     -  assignment
%     alg   -  'psm'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  A. Egozi
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 12-22-2014

    [N1, N2] = size(X0);
    ro = 1e-6;
    iter = 1000;
    p = ones(N1 * N2, 1) / N2;
%    AA = full(K);
    AA = K;
    
    for i = 1 : iter
        p_ = p;
        q = AA * p;
        if rem(i, 2) == 1
            p = normalize1(q, N1, N2);
        else
            p = normalize2(q, N1, N2);
        end
        p_ = p_ + eps;
        
        % modified by Tao Wang, to reduce memory cost
      %  AA = AA .* ((p./ p_) * ones(1, N1 * N2));
        if issparse(AA)
            scale = p ./ p_;
            [row, col, v] = find(AA);
            for k = 1:size(row)
                AA(row(k), col(k)) = AA(row(k), col(k)) * scale(row(k));
            end
        else
            AA = AA .* ((p./ p_) * ones(1, N1 * N2));
        end
        
        if sqrt(sum((p - p_).^2))/(N1 * N2) < ro
            break;   
        end 
    end
    X = reshape(p, N2, N1);     
    X = X';
end

function p = normalize1(q, N1, N2)
% normalized by row

% modified by Tao Wang, to reduce memory cost

% norm_q = kron(eye(N1), ones(N2, N2)) * q;
% p = q ./ norm_q;

q = reshape(q, N2, N1)';
norm_q = sum(q, 2) + eps;    % sum by row
for i = 1:N1
    q(i,:) = q(i,:) / norm_q(i);
end
p = reshape(q', N1*N2, 1);
end

function p = normalize2(q, N1, N2)
% normalized by row

% modified by Tao Wang, to reduce memory cost
% norm_q_2 = kron(ones(N1, N1), eye(N2)) * q;
% p = q ./ norm_q_2;

q = reshape(q, N2, N1)';
norm_q = sum(q, 1) + eps;    % sum by row
for i = 1:N2
    q(:,i) = q(:,i) / norm_q(i);
end
p = reshape(q', N1*N2, 1);
end

% function P = permutation(vec, m, n)
%     P = zeros(m, n);
%     for i = 1 : m
%        P(i, vec(i)) = 1; 
%     end   
% end

% function X = gmPosCPsm(K, X0, par)
% % Probabilistic Spectral Matching.
% %
% % Reference
% %   A. Egozi, Y. Keller, and H. Guterman. "A Probabilistic Approach to
% %   Spectral Graph Matching", IEEE Trans. Pattern Anal. Mach. Intell.,
% %   35(1):18-27, 2013.
% %
% % Math
% %   This algorithm is to obtain the optimal X for the following problem
% %     max_x   x' * K * x
% %     s.t.    A * x <= 1
% %
% % Remark
% %   nn = n1 x n2
% %
% % Input
% %   K        -  affinity matrix, nn x nn (sparse)
% %   X0       -  initial assignment, n1 x n2
% %   par      -  parameter
% %
% % Output
% %   X        -  permutation matrix, n1 x n2
% %
% % History
% %   create   -  Tao WANG (twang@bjtu.edu.cn), 12-20-2014
% 
%     % dimension
%     [n1, n2] = size(X0);
%     P0 = ones(n1 * n2, 1) / n2;
% 
%     % function parameter
%     d = 0.0001;
%     
%     K = full(K);
%     Pt = P0;
%   %  Pt = sparse(P0);
%     for t = 1:20
%         Qt = K * Pt;
%         Pnew = Normalize(Qt, n1, n2); 
%         
%         for i = 1:n1*n2
%             for k = 1:n1*n2
%                 if K(i, k) > 0.0 && Pt(i) > 0.0
%                     K(i, k) = K(i, k) * Pnew(i) / Pt(i);
%                 end
%             end
%         end
%        
%         deta = Pnew - Pt;
%       %  deta = full(Pnew - Pt);
%         err = norm(deta) / (n1 * n2);
%         if (err < d) 
%             break;
%         end
%         Pt = Pnew;
%     end
% 
%     X = reshape(Pnew, [n1 n2]);
% end
% 
% function P = Normalize(P0, n1, n2)
%     X = reshape(P0, [n1 n2]);
%     S = sum(X, 2);
%     for i = 1:n1
%         for k = 1:n2
%             X(i,k) = X(i,k) / S(i);
%         end
%     end
%     P = reshape(X, n1 * n2, 1);
% end