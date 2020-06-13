function X = GNCCP_K(K, group1, group2, par)
% graph matching by GNCCP
%
% Reference
%[1] Zhi-Yong Liu and Hong Qiao, "GNCCP - Graduated NonCovexity and Concavity
%Procedre", IEEE Trans. PAMI, DOI: 10.1109/TPAMI.2013.223, 2014
%[2] Zhi-Yong Liu, H. Qiao, X. Yang and Steven C.H. Hoi, "Graph Maching by
%Simplified Convex-Concave Relaxation Procedure", IJCV, DOI: 10.1007/s11263-014-0707-7, 2014
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
%   X: steady state distribution of GNCCP
%
% History
%   create  -  Zhiyong Liu
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 04-11-2015
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 05-05-2015
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 07-25-2015
%           -  code optimization for rapid computation

%% generate match info
[mid1, gid1] = find(group1);
[mid2, gid2] = find(group2);
nGroups1 = max(gid1);
nGroups2 = max(gid2);
nMatches = size(K, 1);
matchInfo = zeros(nMatches, 2);
matchInfo(mid1, 1) = gid1;
matchInfo(mid2, 2) = gid2;

%% parameter
eta = 0.001; %for gradient value
deta = ps(par, 'deta', 0.01);
nItMa = ps(par, 'nItMa', 100);


%% initialize
gamma = 1;
K = -K;
x = ones(nMatches, 1);
if issparse(K)
    x = sparse(x);
end
x = normalize(x, group1);

KD = K'+K;  % store for rapid computation

while gamma > -1
    
    % Frank-wolfe
    x = FW(x, K, KD, matchInfo, nGroups1, nGroups2, gamma, nItMa, eta);
   
    if max(nGroups1, nGroups2) - x'*x < eta
        break;
    end
    
    gamma = gamma - deta;
%     if mod(round(gamma * 1000), 100) == 0
%         disp(['gamma=', num2str(gamma)]);
%     end
end %end of gamma

X = x;

end

function x = FW(x, K, KD, matchInfo, nGroups1, nGroups2, gamma, nItMa, eta)
    nMatches = size(K, 1);
    for ite = 1:nItMa
        if gamma > 0
        %    g = 2 * gamma * x + (1-gamma) * (K' + K) * x;
            g = 2 * gamma * x + (1-gamma) * (KD * x);
        else
        %    g = 2 * gamma * x + (1+gamma) * (K' + K) * x;
            g = 2 * gamma * x + (1+gamma) * (KD * x);
        end
        
        if (g' * g) < eta
            break;
        end
        
       %% subproblem 
       %% y = reshape(KM(-reshape(g,n1,n2)),n1 * n2,1);
        maxv = max(g); 
        G = zeros(max(nGroups1, nGroups2));
        G(:,:) = maxv + 1;
        for mid = 1:nMatches
            G(matchInfo(mid,1), matchInfo(mid, 2)) = g(mid);
        end
        G = KM(-G);
        y = zeros(nMatches, 1);
        for mid = 1:nMatches
            y(mid) = G(matchInfo(mid,1), matchInfo(mid, 2));
        end
        
        if issparse(K)
            y = sparse(y);
        end
        errorterm = sum(sum(g.*(y-x)));
        if(errorterm > -eta) 
            break; 
        end
        
        %line search, divided method 1
%         alpha = linesearch(y,x,K,gamma);
%         x = x + alpha * (y - x);
        x = linesearch2(y, x, K, gamma, eta);
 
        Fnew = fl(x, K, gamma);
        if gamma > 0
        %    g = 2 * gamma * x + (1-gamma) * (K' + K) * x;
            g = 2 * gamma * x + (1-gamma) * (KD * x);
        else
        %    g = 2 * gamma * x + (1+gamma) * (K' + K) * x;
            g = 2 * gamma * x + (1+gamma) * (KD * x);
        end
        
        tmp = sum(sum(g.*(x - y)));
        if tmp < eta  * abs(Fnew - tmp);
            break;
        end      
    end %end of frank-wolfe

end

function f = fl(x, K, gamma)
    if gamma > 0
        f = gamma * (x'*x) + (1-gamma) * (x' * K * x );
    else
        f = gamma * (x'*x) + (1+gamma) * (x' * K * x);
    end
end

function x = linesearch2(y, x, K, gamma, eta)
    Pright = y;
    Pleft = x;
    deltaX = (Pright - Pleft);
    for i = 1:10
        Pnew = Pleft + 0.5 * deltaX;
        Pnew2 = Pleft + (0.5 + eta) * deltaX;
        F0 = fl(Pnew, K, gamma);
        F1 = fl(Pnew2, K, gamma);

        if F0 < F1
           Pright = Pnew;
        else
           Pleft = Pnew;
        end
        deltaX = Pright - Pleft;
    end
    x = Pnew;
    %%
end


function alpha = linesearch(y,x,K,gamma)
z = y-x;
t1 = z'*K*z;
t2 = z'*z;
t3 = x'*K*z;
t4 = z'*K*x;
t5 = z'*x;
t6 = x'*z;
if gamma > 0
    a = (1-gamma) * t1 + gamma * t2;
    b = (1-gamma) * (t3 + t4) + gamma * (t5 + t6);
else
    a = (1+gamma) * t1 + gamma * t2;
    b = (1+gamma) * (t3 + t4) + gamma * (t5 + t6);
end

if a > 0
    alpha = - 0.5 * b/a;
    if alpha > 1
        alpha = 1;
    end
    if alpha < 0
        alpha = 0;
    end
else
    if a+b >= 0
        alpha = 0;
    else
        alpha = 1;
    end
end

end
