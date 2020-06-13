function X = BPF_G(K, group1, group2, paraStr)
% extension based on GNCCP
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
%   create  -  Tao WANG (twang@bjtu.edu.cn), 05-21-2015
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 07-25-2015
%           -  code optimization for rapid computation
%   modify  -  Tao WANG (twang@bjtu.edu.cn), 10-20-2016
%           -  combine APE for improvement

% store for rapid computation
global Knew KD groupMin groupMax;

%% parameter
para.nItMa  = ps(paraStr, 'nItMa', 100);
para.deta = ps(paraStr, 'deta', 0.01);
para.nHist = ps(paraStr, 'nHist', 10);
para.eta = ps(paraStr, 'eta', 1e-3);        % for FW algorithm thresold
para.tol = ps(paraStr, 'tol', 1e-3);        % for arclength function tol
para.theta = ps(paraStr, 'theta', 1e-3);    % for adaptive path following
para.rho = ps(paraStr, 'rho', 2);        % for adaptive path following

%% generate match group info
group1 = logical((group1));
group2 = logical((group2));
if size(group1, 2) <= size(group2, 2)
    groupMin = group1;
    groupMax = group2;
else
    groupMin = group2;
    groupMax = group1;
end
nMatches = size(K,1);

%% initialize
Knew = -K;
KD = Knew' + Knew;
gamma = 1;
x = ones(nMatches, 1);
if issparse(K)
    x = sparse(x);
end
x = normalize(x, group1);

% original path following
[path, steps, x] = PathFollowing_Adaptive(x, gamma, group1, group2, para);
x = greedyMapping(x, group1, group2);
optX   = x;
maxObj = x' * K * x;
path_count  = 0;
sp_count    = 0;
fprintf('\n Objective of original path: %.4f', full(maxObj));

gamma_arc = 1;

while 1
    gamma = gamma_arc;
    
    if issparse(K)
        u_path = sparse(nMatches + size(groupMin, 2) + size(groupMax, 2), size(path, 2));
    else
        u_path = zeros(nMatches + size(groupMin, 2) + size(groupMax, 2), size(path, 2));
    end
    s_Jaco = zeros(1, size(path, 2));
    [alpha, beta, mu, ~] = computeCacheMatrix(nMatches, groupMin, groupMax);

    
    bNewPath = 0;
    for i = 1 : size(path, 2) - 1
        x = path(:, i);
        
        [flag, alpha, beta, mu, jaco] = KKT_Multiplier_iter(x, gamma, alpha, beta, mu, para);
        s_Jaco(i) = jaco;
        u_path(1:nMatches, i) = x;
        u_path(nMatches+1:end, i) = [alpha; beta];
        if ~flag
            break;
        end
       
        i_jaco = i - 1;
        if isSingularPoint(s_Jaco, i_jaco)
            sp_count = sp_count + 1;

            % tangent vector estimated by previous points
            k = max(1, i_jaco - para.nHist + 1);
            steps_store = steps(k:i_jaco);
            u_gra = EstimateGradient(u_path(:, k : i_jaco), steps_store);
        %    u_gra = EstimateGradient(u_path(:, k : i_jaco));
            
            % branching a new path
        %    gamma_jaco = gamma + para.deta;
            gamma_jaco = gamma + steps(i_jaco) * para.deta;
            [branched, gamma_arc, x_arc] = branch_by_arc(gamma_jaco, i_jaco, steps(i_jaco), u_path, u_gra, para);
            if ~branched
                continue;
            end
           
            % compute the path using branched point    
            [path_arc, steps_arc, x] = PathFollowing_Adaptive(x_arc, gamma_arc, group1, group2, para);
            x = greedyMapping(x, group1, group2);
            obj = x' * K * x;
 
            if obj > maxObj
                maxObj   = obj;
                optX     = x;
                optPath  = path_arc;
                
                bNewPath = 1;
                path_count = path_count + 1;
                fprintf('\n\n  ======= gamma = %.4f, branched a better path, objectives ========= %.4f', gamma, full(maxObj));
                break;
            end
        end    

        gamma = gamma - steps(i) * para.deta;
     %   gamma = gamma - para.deta;    
    end
    
    if ~bNewPath
        break;
    end

    path = optPath;
    steps = steps_arc;
end    

X = optX;

end

function [branched, gamma_arc, x_arc] = branch_by_arc(gamma, ind_jaco, step_jaco, u_path, u_gra, para)
global groupMin groupMax;

nMatches = size(groupMin, 1);
%step_len = 1;
step_len = step_jaco;
branched = 0;
gamma_arc = gamma - step_len * para.deta;
x_arc = [];

if gamma_arc > -1 
    u0 = u_path(:, ind_jaco);
    u_tan = step_len * u_gra;
    [flag, u_arc] = arclength(gamma_arc, [groupMin, groupMax], u0, u_tan, para);
 
    if ~flag
        fprintf('\n\t gamma = %.4f, fail to compute arc length !!!!!!!!!', gamma_arc);
        return;
    end
       
    x_arc = u_arc(1:nMatches);
    x_arc(x_arc<0) = 0;
    x_arc(x_arc>1) = 1;
    x_arc = normalize(x_arc, groupMin);

    % FW
    x_arc = FW(x_arc, groupMin, groupMax, gamma_arc, para.nItMa, para.eta);
    if (ind_jaco + 1) <= size(u_path, 2)
        x_org = u_path(1:nMatches, ind_jaco + 1);
        if norm(x_arc - x_org) < 1e-3
            fprintf('\n\t gamma = %.4f, invalid singular point!!!!!!!!', gamma_arc);
            return;
        end
    end
    branched = 1;   
end

end

function [flag, alpha, beta, mu, jaco] = KKT_Multiplier_iter(x, gamma0, alpha, beta, mu, para)
%% this function is to compute the KKT multiplier \alpha, \beta, \mu
%% at current solution point x, given initial guess of \alpha, \beta, \mu. 
%% and further compute the sign of the jacobian matrix.

global KD groupMin groupMax;
gamma = gamma0;
if gamma > 0
    graX = 2 * gamma * x + (1 - gamma) * KD * x;
else
    graX = 2 * gamma * x + (1 + gamma) * KD * x;
end
        
optType = 0;

[flag, alpha, beta, mu] = KKT(optType, x, graX, groupMin, groupMax, alpha, beta, mu, para);
[jaco,~] = JDetSign(optType, x, gamma, beta, mu, groupMin, groupMax);
end

function [path, steps, x] = PathFollowing_Adaptive(x, gamma, groupMin, groupMax, par)
global Knew;

% initial guess and solutions
nGroupsMax = size(groupMax, 2);
nMatches = size(groupMax, 1);
nHist = par.nHist;

path = [];
steps = [];
count = 0;
sLen = 1;

% storage of history temporal solutions
if issparse(x)
    x_store = sparse(nMatches, nHist);
else
    x_store = zeros(nMatches, nHist);
end
sLen_store = zeros(1, nHist);


while gamma > -1
    x = FW(x, groupMin, groupMax, gamma, par.nItMa, par.eta);
    
    path = [path, x];
    
    %% check for convergency of gamma
    if nGroupsMax - x'*x < par.eta
        break;
    end
   
    % adjust step length
    if count >= nHist
        sLen = AdaptiveStep(x, x_est, sLen, par.theta, par.rho);
    end
    
    steps = [steps, sLen];
    gamma = gamma - sLen * par.deta;

    % store history solutions
    x_store(:, 1:nHist-1) = x_store(:, 2:nHist);
    x_store(:, nHist) = x;
    
    % store history step length
    sLen_store(1:nHist-1) = sLen_store(2:nHist);
    sLen_store(nHist) = sLen;

    % estimate next solution x as the start point of next iteraton
    count = count + 1;
    if count >= nHist
       x_est = EstimateByTangent(x_store, sLen_store, 1);
       
       x_est(x_est<0)=0;
       x_est(x_est>1)=1;
       x_est = normalize(x_est, groupMin);
        
      if fl(x_est, Knew, gamma) < fl(x, Knew, gamma)
            x = x_est;
      end
    end

end


end

function [flag, alpha, beta, mu] = KKT(optType, x, graX, G1, G2, alpha0, beta0, mu0, para)

nn = size(x,1);
n1 = size(G1, 2);
n2 = size(G2, 2);

flag = 1;

%% number of parameters in \mu to be computed 
alpha = alpha0;
if issparse(mu0)
    mu = sparse(nn, 1);
else
    mu    = ones(nn, 1);
end
mu(x>0) = 0;

if n1 == n2
    beta = beta0;
else
    if issparse(beta0)
        beta = sparse(n2, 1);
    else
        beta  = ones(n2, 1);
    end
    sum_g2 = x'*G2;
    beta(sum_g2 < 1) = 0;
end

%if sum(mu) + sum(beta) > 0
%if (n1 == n2 && sum(mu) > 0) || (n1 < n2 && sum(mu) + sum(beta) > 0)
    dmu = diag(-mu);
    dmu = dmu(:,any(dmu));
    if (n1 < n2)
        dbeta = G2(:, beta>0);
        if optType
            A = [G1, -dbeta, -dmu];
        else
            A = [G1, dbeta, dmu];
        end
    elseif (n1 == n2)
        if optType
            A = [G1, G2, -dmu];
        else
            A = [G1, G2, dmu];
        end
    end
    
    if size(A, 1) < size(A, 2)
        flag = 0;
        return;
    end
    
    if n1 < n2
        y0 = [alpha0; beta0(beta>0); mu0(x<=0)]; 
    elseif n1 == n2
        y0 = [alpha0; beta0; mu0(x<=0)];
    end

    UB = inf(size(A,2), 1);
    LB = -inf(size(A,2), 1);
    if (n1 < n2)
        LB(n1 + 1 : end) = 0;
    elseif (n1 == n2)
        LB(n1 + n2 + 1 : end) = 0;
    end
    
%   maxIter = min(100, ceil(2 * sqrt(size(x,1))));
   maxIter = min(10, ceil(2 * sqrt(size(x,1))));
   options = optimoptions('lsqlin', 'Display', 'off','TolFun', para.tol, 'MaxIter', maxIter);
   y = lsqlin(A, -graX, [], [], [], [], LB, UB, y0, options);
    
    if any(isnan(y)) || any(isinf(y))
        fprintf('\n invalid solution of KKT equation system!!!!!!!!!!!!!!!');
        flag = 0;
        return;
    end
    
    alpha  = y(1:n1);
    if n1 < n2
        beta(beta>0)= y(n1+1 : n1+size(dbeta,2));
        mu(mu>0)    = y(n1+size(dbeta,2)+1 : end);
        
        assert(all(beta >= 0));
    elseif n1 == n2
        beta = y(n1 + 1 : n1 + n2);
        mu(mu>0) = y(n1 + n2 + 1 : end);
    end

    assert(all(mu >= 0));
%end

end

function [flag, u_arc] = arclength(gamma_arc, G, u0, u_tan, para)
global KD;
%global groupMin groupMax;
nn = size(KD, 1);
nG = size(G, 2);
%n1 = size(groupMin, 2);
%n2 = size(groupMax, 2);

if issparse(KD)
    I = speye(nn);
    ZnG = sparse(nG, nG);
    Znn = sparse(nn, 1);
    OnG = sparse(ones(nG, 1));
else
    I = eye(nn);
    ZnG = zeros(nG, nG);
    Znn = zeros(nn, 1);
    OnG = ones(nG, 1);
end

if gamma_arc > 0
    coe = 2 * gamma_arc * I + (1 - gamma_arc) * KD;
else
    coe = 2 * gamma_arc * I + (1 + gamma_arc) * KD;
end

A = [coe, G;
     G', ZnG;
    u_tan'];
b = [Znn; OnG; u_tan' * (u_tan + u0)];
    
UB = inf(size(A,2), 1);
LB = -inf(size(A,2), 1);
LB(1 : nn) = 0;

if issparse(u0)
    UB = sparse(UB);
    LB = sparse(LB);
end

u0 = u0 + u_tan;
maxIter = min(10, ceil(2 * sqrt(nn)));
options = optimoptions('lsqlin', 'Display', 'off','TolFun', para.tol, 'MaxIter', maxIter);
u_arc = lsqlin(A, b, [], [], [], [], LB, UB, u0, options);

if any(isnan(u_arc)) || any(isinf(u_arc))
    u_arc = u0;
    flag = 0;
else
    flag = 1;
end

if issparse(u0)
    u_arc = sparse(u_arc);
end

end

function [s, l] = JDetSign(optType, x, gamma, beta, mu, G1, G2)
global KD;

n1 = size(G1, 2);
n2 = size(G2, 2);
nn = size(x, 1);

if issparse(x)
    I = speye(nn);
    Zn1n1 = sparse(n1, n1);
    Zn1n2 = sparse(n1, n2);
    Zn1nn = sparse(n1, nn);
    Zn2n1 = sparse(n2, n1);
    Zn2n2 = sparse(n2, n2);
    Zn2nn = sparse(n2, nn);
    Znnn1 = sparse(nn, n1);
    Znnn2 = sparse(nn, n2);    
else
    I = eye(nn);
    Zn1n1 = zeros(n1, n1);
    Zn1n2 = zeros(n1, n2);
    Zn1nn = zeros(n1, nn);
    Zn2n1 = zeros(n2, n1);
    Zn2n2 = zeros(n2, n2);
    Zn2nn = zeros(n2, nn);
    Znnn1 = zeros(nn, n1);
    Znnn2 = zeros(nn, n2);    
end

if gamma > 0
    coe = 2 * gamma * I + (1 - gamma) * KD ;
else
    coe = 2 * gamma * I + (1 + gamma) * KD;
end

if (n1 < n2)
    if issparse(x)
        B = beta * sparse(ones(1, nn));
    else
        B = beta * ones(1, nn);
    end
    C = x' * G2 - 1;
    C = diag(C);
    if optType
        J = [coe, G1, -G2, I; ...
            G1', Zn1n1, Zn1n2, Zn1nn; ...
            G2' .* B, Zn2n1, C, Zn2nn; ...
            diag(mu), Znnn1, Znnn2, diag(x)];
    else
        J = [coe, G1, G2, -I; ...
            G1', Zn1n1, Zn1n2, Zn1nn; ...
            G2' .* B, Zn2n1, C, Zn2nn; ...
            diag(mu), Znnn1, Znnn2, diag(x)];
    end
elseif (n1 == n2)
    if optType
        J = [coe, G1, G2, I; ...
            G1', Zn1n1, Zn1n2, Zn1nn; ...
            G2', Zn2n1, Zn2n2, Zn2nn; ...
            diag(mu), Znnn1, Znnn2, diag(x)];
    else
        J = [coe, G1, G2, -I; ...
            G1', Zn1n1, Zn1n2, Zn1nn; ...
            G2', Zn2n1, Zn2n2, Zn2nn; ...
            diag(mu), Znnn1, Znnn2, diag(x)];
    end    
end

%[~, U, P] = lu(J);
[~, U, P, Q] = lu(J);

du = diag(U);
%dl = diag(L);
%s = prod(sign(du)) * prod(sign(P)) * prod(sign(Q));
s = det(Q) * det(P) * prod(sign(du));
%s = prod(sign(du));
%l = sum(log(abs(du)));
l = 0;

end

function [alpha, beta, mu, gra] = computeCacheMatrix(nMatches, groupMin, groupMax)
n1 = size(groupMin, 2);
n2 = size(groupMax, 2);
% alpha = zeros(n1, 1);
% beta  = zeros(n2, 1);
% mu    = zeros(nMatches, 1);
% gra   = zeros(2 * nMatches + n1 + n2, 1);

alpha = sparse(n1, 1);
beta  = sparse(n2, 1);
mu    = sparse(nMatches, 1);
gra   = sparse(2 * nMatches + n1 + n2, 1);


end


function flag = isSingularPoint(s_Jaco, i)
%% we regard there is a singular point if two adjcent point have Jacobian with different sign
if i <= 1 || i >= size(s_Jaco, 2)
    flag = 0;
    return;
end

if s_Jaco(i) * s_Jaco(i + 1) <= 0
    flag = 1;
else
    flag = 0;
end

end

function gra_mean = EstimateGradient(x_store, steps_store)
%% estimate the gradient at current point according to previous K solutions
% x_store:  [N * K] matrix
if issparse(x_store)
    gra_mean = sparse(size(x_store, 1), 1);
else
    gra_mean = zeros(size(x_store, 1), 1);
end
k = size(x_store, 2);
for i = 1:k-1
    gra = (x_store(:,i+1) - x_store(:,i)) / steps_store(i);
    gra_mean = gra_mean + i * gra;
end
gra_mean = gra_mean / sum(1:k-1);
end


function x = FW(x, group1, group2, gamma, nItMa, eta)
global Knew KD;

nMatches = size(KD, 1);
nGroups1 = size(group1, 2);
nGroups2 = size(group2, 2);

for ite = 1:nItMa
    if gamma > 0
        g = 2 * gamma * x + (1-gamma) * (KD * x);
    else
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
        gid1 = group1(mid,:);
        gid2 = group2(mid,:);
        G(gid1, gid2) = g(mid);
    end
    G = KM(-G);
    y = zeros(nMatches, 1);
    for mid = 1:nMatches
        gid1 = group1(mid,:);
        gid2 = group2(mid,:);
        y(mid) = G(gid1, gid2);
    end
        
    if issparse(Knew)
        y = sparse(y);
    end
        
    errorterm = sum(sum(g.*(y-x)));
    if(errorterm > -eta) 
        break; 
    end
    
    
    %% line search, divided method 2
    x = linesearch2(y, x, Knew, gamma, eta);
 
    if gamma > 0
        g = 2 * gamma * x + (1-gamma) * (KD * x);
    else
        g = 2 * gamma * x + (1+gamma) * (KD * x);
    end
        
    %% stop condition
    Fnew = fl(x, Knew, gamma);
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
