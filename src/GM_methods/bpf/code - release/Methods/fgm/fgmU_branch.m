function asg = fgmU_branch(KP, KQ, Ct, gphs, asgT, par)
% Factorized graph matching.
%
% Remark
%   The edge is undirected and the edge feature is symmetric.
%
% Reference
%   F. Zhou and F. De la Torre, "Factorized Graph Matching", in CVPR, 2012.
%
% Input
%   KP        -  node affinity matrix, n1 x n2
%   KQ        -  edge affinity matrix, m1 x m2
%   Ct        -  correspondence constraint, n1 x n2
%                  Ct_ij = 1: i and j can be matched
%                  Ct_ij = 0: i and j cannot be matched
%   gphs      -  graphs, 1 x 2 (cell)
%     G       -  node-edge adjacency matrix, ni x mi
%     H       -  augment node-edge adjacency matrix, ni x (mi + ni)
%   asgT      -  ground-truth assignment (can be [])
%   par       -  parameter
%     nAlp    -  #alpha, {100}
%     nItMa   -  #maximum iteration steps for each scale of alpha, {100}
%     nHst    -  #history nodes for modifed FW algorithm, {10}
%     ip      -  flag of using IPFP to improve the algorithm, {'y'} | 'n'
%     thAlp   -  threshold for alpha used for deciding when to start use CCCP, {0}
%     deb     -  flag of debugging, 'y' | {'n'}
%     idxAlp  -  index of alphas that needed to be explored for more details, {[]}
%
% Output
%   asg       -  assignment
%     alg     -  'fgmU'
%     X       -  correspondence matrix, n1 x n2
%     acc     -  accuracy
%     obj     -  objective
%
% History     
%   create    -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

global coeVex coeCav;

% dimension
G1 = gphs{1}.G;
G2 = gphs{2}.G;
H1 = gphs{1}.H;
H2 = gphs{2}.H;
[n1, m1] = size(G1);
[n2, m2] = size(G2);

% make sure the graphs are of the same size
if n1 ~= n2
    [KP, G1, G2, H1, H2, Ct] = makeEq(KP, G1, G2, H1, H2, Ct);
end
n = max(n1, n2);

% figure for debugging
% if isDeb
%     rows = 1; cols = 2;
%     Ax = iniAx(10, rows, cols, [250 * rows, 250 * cols]);
%     ha = [];
% end

% factorize
fact(KP, KQ, G1, G2, H1, H2);

% initialization
X0 = gmIniUnif(Ct, st('nor', 'doub'));

coeVex = coeConvex(H1, H2, X0);
coeCav = coeConcav(G1, G2, X0);


initSolutions = reshape(X0', n * n ,1);
initGuesses   = 1;

% allocate space for saving variables
%[objs, objRs] = zeross(1, nAlp);
%Xs = cell(1, nAlp);

% path-following
%prCIn('path', nAlp, .1);

%prCOut(nAlp);

optX   = zeros(n, n);
maxObj = 0;
count  = 0;

while size(initSolutions, 2) > 0 && size(initGuesses, 2) > 0
    fprintf('\n ===== Path following from a new guess at iAlp = %d ======', initGuesses(1));
    
    [X, initSolutions, initGuesses] = PathFollowing(n, initSolutions, initGuesses, par);
    
    % save best solution X
    X = gmPosDHun(X);
    [~, obj] = evalObj(X, 1);
    if obj > maxObj
        optX = X;
        maxObj = obj;
    end
    
    count = count + 1;
end

X = optX;

fprintf('\n ** Totaolly %d branches are searched! ** \n', count);

% post-processing (check the integrity)
XC = X;
X = gmPosDHun(X);
if ~equal('XC', XC, X, 'pr', 'n')
    pr('non-discrete');
end

% objective
[~, obj] = evalObj(X, 1);

% re-size to the original size
X = X(1 : n1, 1 : n2);

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% store
asg.alg = 'fgmU';
asg.X = X;
asg.acc = acc;
asg.obj = obj;
%asg.objs = objs;
%asg.objRs = objRs;

prOut;


function [X, initSolutions, initGuesses] = PathFollowing(n, initSolutions, initGuesses, par)
% function parameter
nAlp = ps(par, 'nAlp', 100);
nItMa = ps(par, 'nItMa', 100);
nHst = ps(par, 'nHst', 10);
thAlp = ps(par, 'thAlp', 0);
isIp = psY(par, 'ip', 'n');
isDeb = psY(par, 'deb', 'n');
idxAlp = ps(par, 'idxAlp', []);
%steps = ps(par, 'steps', 1);
prIn('fgmU', 'nAlp %d, nItMa %d, nHst %d, isIp %d', ...
     nAlp, nItMa, nHst, isIp);

% weight
alps = linspace(0, 1, nAlp);

% allocate space for saving variables
[objs, objRs] = zeross(1, nAlp);
%Xs = cell(1, nAlp);


% initial guess and solutions
X0 = reshape(initSolutions(:,1), n, n)';
initAlp = initGuesses(1);

initSolutions = initSolutions(:, 2:end);
initGuesses   = initGuesses(2:end);


% record history solutions and gradients
nHist = ps(par, 'nHist', 5);
x_store = zeros(n * n, nHist);
y_store = zeros(n * n + 2 * n, nHist);
g_store = zeros(n * n, nHist);
g_err   = zeros(1, nAlp - initAlp + 1);

[B, beta, mu, ~] = computeCacheMatrix(n);

kkt_flag = 1;
for iAlp = initAlp : nAlp
    if mod(iAlp, 100) == 0
        fprintf('\n FGM_U: iAlp = %d', iAlp);
    end
%    prC(iAlp);

    % scale of alpha
    alp = alps(iAlp);

    % FW
    [X,~,~,~] = mfw(X0, alp, nItMa, nHst, isDeb);
    
    % CCCP
    if alp < thAlp
        X = cfw(X0, alp, 10, 20, nHst, isDeb);
    end

    % IPFP
    if isIp
        % objective
        [objs(iAlp), objRs(iAlp)] = evalObj(X, alp);

        % using IPFP to refine the result
        if iAlp > 1 && objRs(iAlp) < objRs(iAlp - 1)
            XIp = ipfpSStep(H1, H2, X0);
            X = XIp;
        end
    end

    if kkt_flag
        x = reshape(X', n * n, 1);

        GrVex = gradVex(X);
        GrCav = gradCav(X);
        graX   = (1 - alp) * GrVex + alp * GrCav;
        graX   = reshape(graX', n * n, 1);
        [flag, beta, mu] = KKT(x, graX, B, beta, mu);

        if ~flag
            fprintf('\n iAlp = %d, failure to compute KKT system, close flag!', iAlp);
            kkt_flag = 0;
        end
    end
    
    if kkt_flag
        % record history solutioins and gradients
        x_store(:, 1 : nHist - 1) = x_store(:, 2 : nHist);
        x_store(:, nHist) = x;

        if iAlp > initAlp
            g = (x_store(:, nHist) - x_store(:, nHist - 1)) * nAlp;
            g_store(:, 1 : nHist - 1) = g_store(:, 2 : nHist);
            g_store(:, nHist) = g;

            g_err(iAlp - initAlp + 1) = norm(g - g_store(:, nHist - 1));
        end

        if iAlp - initAlp > nHist && isSingularPoint(g_err(1:iAlp - initAlp + 1), nHist)
                % add a new guess using the estimated solution
            steps = 1;
            branched = 0;
            while (iAlp - 1 + steps) <= nAlp && ~branched
                y0 = y_store(:, nHist);
                [flag, y_arc] = arclength(n, y0, iAlp - 1, steps, B, y_store);
                if ~flag
                    fprintf('\niAlp = %d, fail to compute arc length !!!', iAlp);
                else
                    x_arc = full(y_arc(1:n*n));
                    x_arc(x_arc<0) = 0;
                    x_arc(x_arc>1) = 1;
                    X_arc = reshape(x_arc, n, n)';
                    X_arc = bistocNormalize_slack(X_arc, 1e-3);

                    % FW
                    alp_arc = alps(iAlp - 1 + steps);
                    [X_arc,~,~,~] = mfw(X_arc, alp_arc, nItMa, nHst, isDeb);

                    if norm(X_arc - X) > 1e-4
                        initGuesses   = [initGuesses, iAlp - 1 + steps];
                        initSolutions = [initSolutions, reshape(X_arc', n * n, 1)];
                        fprintf('\n iAlp = %d, singular point detected, add a new guess!!!!!!', iAlp);
                        branched = 1;
                    end
                end

                steps = steps * 2;
            end
        end    

        y = [x; beta];
        y_store(:, 1 : nHist - 1) = y_store(:, 2 : nHist);
        y_store(:, nHist) = y;
    end
    
    % set start point of next iteration
    X0 = X;
end



function [flag, y] = arclength(n, y0, alp0, steps, G, y_store)
global coeVex coeCav;
global GKGP;

%f'(x) = coe * x + extra
alp_arc = alp0 + steps;
coe     = (1 - alp_arc) * coeVex + alp_arc * coeCav;
extra   = alp_arc * reshape(GKGP', n * n, 1);

nn = size(coe,1);
nG = size(G, 2);

y_tan = steps * EstimateGradient(y_store);
A = [coe, G;
     G', sparse(nG, nG);
    sparse(y_tan')];
b = [-extra; ones(nG, 1); y_tan' * (y_tan + y0)];
b = sparse(b);
    
UB = inf(size(A,2), 1);
LB = -inf(size(A,2), 1);
LB(1 : nn) = 0;
options = optimoptions('lsqlin', 'Display', 'off');
y0 = y0 + y_tan;
y = lsqlin(A, b, [], [], [], [], LB, UB, y0, options);
%y = A \ b;

fprintf('\n y_tan = %.4f, y_err = %.4f, kkt_err = %.4f', ...
     norm(y_tan), norm(y-y0), norm(A*y-b));

if any(isnan(y)) || any(isinf(y))
    y = y0;
    flag = 0;
else
    flag = 1;
end

%x = y(1:nn);
%assert(all(x >=0));

function gra = EstimateGradient(x_store)
%% estimate next solution x according to previous K solutions
% x_store:  [N * K] matrix
gra = zeros(size(x_store, 1), 1);
k = size(x_store, 2);
for i = 1:k-1
    shift = (x_store(:,i+1) - x_store(:,i));
    gra = gra + i * shift;
end
gra = gra / sum(1:k-1);


function flag = isSingularPoint(g_err, nHist)
%% we regard there is a singular point if the gradient error increases suddenly
flag = 0;
k = size(g_err, 2);
if k > nHist
    err_est = sum([1:nHist] .* g_err(k - nHist: k - 1)) / sum(1:nHist);
    if g_err(k) > 5 * err_est
        flag = 1;
    end
end


function x = PathEstimate(x_store, steps)
%% estimate next solution x according to previous K solutions
% x_store:  [N * K] matrix
% scale_store: K vector

x_shift = zeros(size(x_store));

%% estimate shift of next solution
k = size(x_store, 2);
for i = 1:k-1
    x_shift(:,i) = (x_store(:,i+1) - x_store(:,i));
    x_shift(:,k) = x_shift(:,k) + i * x_shift(:,i);
end
x_shift(:,k) = x_shift(:,k) / sum(1:k-1);

x = x_store(:,k) + x_shift(:,k) * steps;

x(x<0)=0;
x(x>1)=1;


function [flag, beta, mu] = KKT(x, graX, G, beta0, mu0)

nn = size(x,1);

flag = 1;

%% number of parameters in \mu to be computed 
beta = beta0;
mu    = sparse(ones(nn, 1));
mu(x>0) = 0;

nG = size(G, 2);

%if sum(mu) > 0 
    dmu = diag(-mu);
    dmu = dmu(:,any(dmu));

    optType = 1;    % maximize objective function
    if optType
        A = [G, -dmu];
    else
        A = [G, dmu];
    end
    
    if size(A, 1) < size(A, 2)
        flag = 0;
        return;
    end

    UB = inf(size(A,2), 1);
    LB = -inf(size(A,2), 1);
    LB(nG : end) = 0;
    
   y0 = sparse([beta0; mu0(x<=0)]);
   maxIter = min(100, ceil(2 * sqrt(size(x,1))));
   options = optimoptions('lsqlin', 'Display', 'off','TolFun', 1e-5, 'MaxIter', maxIter);
   y = lsqlin(A, -graX, [], [], [], [], LB, UB, y0, options);
%     options = optimoptions('lsqlin', 'Display', 'off');
%     y = lsqlin(A, -graX, [], [], [], [], LB, UB, [], options);
    
    if any(isnan(y)) || any(isinf(y))
        flag = 0;
        return;
    end
    
    beta  = y(1:nG);
    mu(mu>0) = y(nG + 1 : end);

    assert(all(mu >= 0));
%end

function [B, alpha0, mu0, gra0] = computeCacheMatrix(n)
Grp1 = sparse(n * n, n);
Grp2 = sparse(n * n, n);
for i = 1:n
    g1 = sparse(n, n);
    g2 = sparse(n, n);
    g1(i,:) = 1;
    g2(:,i) = 1;
    
    Grp1(:,i) = reshape(g1', n * n, 1);
    Grp2(:,i) = reshape(g2', n * n, 1);
end
B = [Grp1, Grp2];
alpha0 = sparse(n + n, 1);
mu0 = sparse(n * n, 1);
gra0 = sparse(2 * n * n + 2 * n, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fact(KP0, KQ0, G1, G2, H1, H2)
% Compute the factorization.

% global variable
global L KP KQ;
global HH1 HH2 UU VV UUHH VVHH HUH HVH GKGP;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;

KP = KP0;
KQ = KQ0;

% L
L = [KQ, -KQ * G2'; -G1 * KQ, G1 * KQ * G2' + KP];

% SVD
[U, S, V] = svd(L);
s = diag(S);
idx = 1 : length(s);
%idx = find(s > 0);
pr('svd: %d/%d', length(idx), length(s));
k = length(idx);
U = U(:, idx);
V = V(:, idx);
s = s(idx);

U = multDiag('col', U, real(sqrt(s)));
V = multDiag('col', V, real(sqrt(s)));

% the following decomposition works very badly
% U = eye(size(L, 1));
% V = L;

% auxiliary variables that will be frequently used in the optimization
UU = U * U';
VV = V * V';

HH1 = H1' * H1;
HH2 = H2' * H2;
UUHH = UU .* HH1;
VVHH = VV .* HH2;
HUH = H1 * UUHH * H1';
HVH = H2 * VVHH * H2';
GKGP = -G1 * KQ * G2' + KP;

% index
IndG1 = mat2ind(G1);
IndG2 = mat2ind(G2);
IndG1T = mat2ind(G1');
IndG2T = mat2ind(G2');
IndH1 = mat2ind(H1);
IndH2 = mat2ind(H2);
IndH1T = mat2ind(H1');
IndH2T = mat2ind(H2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, objs, ts, nIt] = fw(X0, alp, nItMa, isDeb)
% Original Frank-Wolfe algorithm.
%
% Input
%   X0     -  initial solution, n1 x n2
%   alp    -  alpha
%   nItMa  -  #maximum iteration number
%   isDeb  -  debug flag, 0 | 1
%
% Output
%   X      -  solution, n1 x n2
%   objs   -  objective, 1 x nItMa
%   ts     -  step size, 1 x nItMa
%   nIt    -  #iteration

[objs, ts] = zeross(1, nItMa);

for nIt = 1 : nItMa
    % gradient
    GrVex = gradVex(X0);
    GrCav = gradCav(X0);
    Gr = (1 - alp) * GrVex + alp * GrCav;

    % hungrian for computing the optimal direction
    YY = gmPosDHun(Gr);
    Y = YY - X0;

    % step size
    [aVex, bVex] = stepSizVex(X0, Y);
    [aCav, bCav] = stepSizCav(X0, Y);
    a = (1 - alp) * aVex + alp * aCav;
    b = (1 - alp) * bVex + alp * bCav;
    t = optStep(a, b);

    % update
    X = X0 + t * Y;

    % debug
    if isDeb
        objs(nIt) = evalObj(X, alp);
        ts(nIt) = t;
    end

    % stop condition
    if norm(X(:) - X0(:)) < eps || t < eps
        break;
    end

    % store
    X0 = X;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, objs, ts, nIt] = mfw(X0, alp, nItMa, nHst, isDeb)
% Modified Frank-wolfe algorithm.
%
% Input
%   X0     -  initial solution, n1 x n2
%   alp    -  alpha
%   nItMa  -  #maximum iteration number
%   nHst   -  #history node
%   isDeb  -  debug flag, 0 | 1
%
% Output
%   X      -  solution, n1 x n2
%   objs   -  objective, 1 x nItMa
%   ts     -  step size, 1 x nItMa

[objs, ts] = zeross(1, nItMa);
Ys = cell(1, nHst);

for nIt = 1 : nItMa

    % gradient
    GrVex = gradVex(X0);
    GrCav = gradCav(X0);
    Gr = (1 - alp) * GrVex + alp * GrCav;
    

    % hungrian for computing the optimal direction
    Y = gmPosDHun(Gr);
    V = Y - X0;
    
    
    % save to history
    pHst = mod(nIt - 1, nHst) + 1;
    Ys{pHst} = Y / nHst;

    % alternative direction
    if nIt >= nHst
        W = -X0;
        for iHst = 1 : nHst
            W = W + Ys{iHst};
        end

        vV = multTr(Gr .* V) / norm(V, 'fro');
        vW = multTr(Gr .* W) / norm(W, 'fro');
        if vW > vV
            V = W;
            Ys{pHst} = Y / nHst;
        end
    end

    
   
    % step size
    [aVex, bVex] = stepSizVex(X0, V);
    [aCav, bCav] = stepSizCav(X0, V);
    a = (1 - alp) * aVex + alp * aCav;
    b = (1 - alp) * bVex + alp * bCav;
    t = optStep(a, b);

    % update
    X = X0 + t * V;
    
    % debug
    if isDeb
        objs(nIt) = evalObj(X, alp);
        ts(nIt) = t;
    end

   
    % stop condition
    % modified by Tao Wang(twang@bjtu.edu.cn), 07/29/2015
%    if norm(X(:) - X0(:)) < eps || t < eps
    eta = 0.001;
%    if norm(X(:) - X0(:)) < eta || t < eta
    if norm(X(:) - X0(:)) * norm(X(:) - X0(:)) < eta || t < eta
        break;
    end

    
    % store
    X0 = X;
end

%pr('mfw, nIt %d', nIt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, objs, objRs] = cfw(X0, alp, nItOutMa, nItInMa, nHst, isDeb)
% CCCP + Frank-wolfe algorithm.
%
% Input
%   X0        -  initial solution, n1 x n2
%   alp       -  alpha
%   nItOutMa  -  #maximum iteration number
%   nItInMa   -  #maximum iteration number
%   nHst      -  #history node
%   isDeb     -  debug flag
%
% Output
%   X         -  solution, n1 x n2
%   objs      -  objective, 1 x nItMa

[objs, objRs] = zeross(1, nItOutMa * nItInMa);
nItDeb = 0;
Ys = cell(1, nHst);

for nItOut = 1 : nItOutMa
    % gradient
    GrCav = gradCav(X0);

    % Frank-Wolfe algorithm for convex optimization
    XX0 = X0;
    for nItIn = 1 : nItInMa

        % gradient
        GrVex = gradVex(XX0);
        Gr = (1 - alp) * GrVex + alp * GrCav;

        % hungrian for computing the optimal direction
        Y = gmPosDHun(Gr);
        V = Y - XX0;

        % save to history
        pHst = mod(nItIn - 1, nHst) + 1;
        Ys{pHst} = Y / nHst;

        % alternative direction
        if nItIn >= nHst
            W = -XX0;
            for iHst = 1 : nHst
                W = W + Ys{iHst};
            end

            vV = multTr(Gr .* V) / norm(V, 'fro');
            vW = multTr(Gr .* W) / norm(W, 'fro');
            if vW > vV
                V = W;
                Ys{pHst} = Y / nHst;
                % fprintf('cfw, hst\n');
            end
        end

        % step size
        [aVex, bVex] = stepSizVex(XX0, V);
        aCav = 0;
        bCav = multTr(GrCav, V);
        a = (1 - alp) * aVex + alp * aCav;
        b = (1 - alp) * bVex + alp * bCav;
        t = optStep(a, b);

        % update
        XX = XX0 + t * V;

        if isDeb
            nItDeb = nItDeb + 1;
            [objs(nItDeb), objRs(nItDeb)] = evalObj(XX, alp);
        end

        % stop condition
        if norm(XX(:) - XX0(:)) < eps || t < eps
            break;
        end

        % store
        XX0 = XX;
    end
    X = XX;

    % stop condition
    if norm(X(:) - X0(:)) < eps
        fprintf('cfw: nItOut %d\n', nItOut);
        break;
    end

    % store
    X0 = X;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [obj, objR, objVex, objCav] = evalObj(X, alp)
% Comupte the objective.
%
% Input
%   X       -  correspondence, n1 x n2
%   alp     -  alpha
%
% Output
%   obj     -  J_alpha
%   objR    -  J_gm
%   objVex  -  J_vex
%   objCav  -  J_cav

% global
global L KQ;
global HH1 HH2 UU VV UUHH VVHH HUH HVH GKGP;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;

% trace(L' * (H1' * X * H2) .^ 2);
tmp1 = multTr(L .* multGXH(IndH1T, X, IndH2) .^ 2);
objR = tmp1;

% trace(U * U' * ((H1' * H1) .* (H1' * X * X' * H1)));
tmp2 = multTr(UU, HH1, multGXH(IndH1T, X * X', IndH1));

% trace(V * V' * ((H2' * H2) .* (H2' * X' * X * H2)));
tmp3 = multTr(VV, HH2, multGXH(IndH2T, X' * X, IndH2));

% convex part
objVex = tmp1 - .5 * tmp2 - .5 * tmp3;

% trace(KQ' * (G1' * X * G2) .^ 2);
tmp1 = multTr(KQ, multGXH(IndG1T, X, IndG2) .^ 2);

% trace((-G1 * KQ * G2' + KP)' * X);
tmp2 = multTr(GKGP, X);

% concave part
objCav = tmp1 + tmp2;

% linear interoplation
obj = (1 - alp) * objVex + alp * objCav;

%%%%%%%%%%%%%%%%%%%%%%%%
function Gr = gradVex(X)
% Compute the gradient of the convex part.

% global
global L; % KQ;
%global HH1 HH2 UU VV UUHH VVHH
global HUH HVH; % GKGP;
%global IndG1;
%global IndG2 IndG1T; 
%global IndG2T;
global IndH1 IndH2 IndH1T IndH2T;
%global GXG ;
global HXH;

%GXG = multGXH(IndG1T, X, IndG2);
HXH = multGXH(IndH1T, X, IndH2);

% 2 * H1 * ((H1' * X * H2) .* L) * H2' - H1 * (HH1 .* UU) * H1' * X - X * H2 * (HH2 .* VV) * H2';
Gr = 2 * multGXH(IndH1, HXH .* L, IndH2T) - HUH * X - X * HVH;


%%%%%%%%%%%%%%%%%%%%%%%%
function Gr = gradCav(X)
% Compute the gradient of the concave part.

% global
%global L;
global KQ;
%global HH1 HH2 UU VV UUHH VVHH HUH HVH;
global GKGP;
global IndG1 IndG2 IndG1T IndG2T; % IndH1 IndH2 IndH1T IndH2T;
global GXG;
%global HXH;

GXG = multGXH(IndG1T, X, IndG2);
%HXH = multGXH(IndH1T, X, IndH2);

Gr = 2 * multGXH(IndG1, GXG .* KQ, IndG2T) + GKGP;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coe = coeConvex(H1, H2, X)
% compute the coefficient matrix of the convex part
% global
global L;
global HUH HVH; 

% convex(X) =  2 * H1 * ((H1' * X * H2) .* L) * H2' - H1 * (HH1 .* UU) * H1' * X - X * H2 * (HH2 .* VV) * H2';
%           = 2 * H1 * ((H1' * X * H2) .* L) * H2' - HUH * X  - X * HVH;
% H1' * X
szX = size(X);
coe = coePrefix(H1', szX);

% (H1' * X * H2)
szX(1) = size(H1',1); 
coe = coeSuffix(H2, szX) * coe;

% (H1' * X * H2) .* L
szX = size(L);
coe = coeProduct(L) * coe;

% H1 * ((H1' * X * H2) .* L)
szX = size(L);
coe = coePrefix(H1, szX) * coe;

% H1 * ((H1' * X * H2) .* L) * H2'
szX(1) = size(H1,1);
coe = coeSuffix(H2', szX) * coe;

coe = 2 * coe - coePrefix(HUH, size(X)) - coeSuffix(HVH, size(X));


function coe = coeConcav(G1, G2, X)
% Compute the coefficient matrix of the concave part.

% global
global KQ;

% concav(X) = 2 * G1 * ((G1' * X * G2) .* KQ) * G2' - G1* KQ * G2' + KP;

% G1' * X
szX = size(X);
coe = coePrefix(G1', szX);

% G1' * X * G2
szX(1) = size(G1', 1);
coe = coeSuffix(G2, szX) * coe;

% (G1' * X * G2) .* KQ
szX = size(KQ);
coe = coeProduct(KQ) * coe;

% G1 * ((G1' * X * G2) .* KQ)
szX = size(KQ);
coe = coePrefix(G1, szX) * coe;

% G1 * ((G1' * X * G2) .* KQ) * G2'
szX(1) = size(G1,1);
coe = coeSuffix(G2', szX) * coe;

coe = 2 * coe;


function coe = coePrefix(A, szX)
% compute the coefficient matrix of vec(A*X)  = coe * vec(X),
% where A is m-by-n matrix and X is n-by-n matrices.
% coe = kron(A,I)
coe = kron(A, speye(szX(2)));

function coe = coeSuffix(A, szX)
% compute the coefficient matrix of vec(X*A)  = coe * vec(X),
% where A is n-by-m matrix and X is n-by-n matrices.
% coe = kron(I, A')
coe = kron(speye(szX(1)), A');

function coe = coeProduct(A)
% compute the coefficient matrix of vec(A.*X) = vec(X.*A) = coe * vex(X),
% where A and X are n-by-n matrices.
% coe = [vec(A)'; zeros(n*n-1,n*n)]
[n1, n2] = size(A);
coe = spdiag(reshape(A',n1*n2,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, b] = stepSizVex(X, Y)
% Obtain the step size for the convex part.

% global
global L; % KQ;
%global HH1 HH2 UU VV;
global UUHH VVHH; % HUH HVH GKGP;
%global IndG1 IndG2 IndG1T IndG2T IndH1;
global IndH2 IndH1T; % IndH2T;
%global GXG;
global HXH;

% auxiliary variables
%GYG = multGXH(IndG1T, Y, IndG2);
H1TY = multGXH(IndH1T, Y, []);
YH2 = multGXH([], Y, IndH2);
H1TX = multGXH(IndH1T, X, []);
XH2 = multGXH([], X, IndH2);
HYH = multGXH([], H1TY, IndH2);

% second-order part
tmp1 = multTr(L .* HYH .^ 2);
tmp2 = multTr(UUHH .* (H1TY * H1TY'));
tmp3 = multTr(VVHH .* (YH2' * YH2));
a = tmp1 - .5 * tmp2 - .5 * tmp3;

% first-order part
tmp1 = multTr(L .* HXH .* HYH);
tmp2 = multTr(UUHH .* (H1TX * H1TY'));
tmp3 = multTr(VVHH .* (XH2' * YH2));
b = 2 * tmp1 - tmp2 - tmp3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, b] = stepSizCav(X, Y)
% Obtain the step size for the concave part.

% global
%global L;
global KQ;
%global HH1 HH2 UU VV UUHH VVHH HUH HVH;
global GKGP;
%global IndG1;
global IndG2 IndG1T; % IndG2T IndH1 IndH2 IndH1T IndH2T;
global GXG; % HXH;

% auxiliary variables
GYG = multGXH(IndG1T, Y, IndG2);
%H1TY = multGXH(IndH1T, Y, []);
%YH2 = multGXH([], Y, IndH2);
%H1TX = multGXH(IndH1T, X, []);
%XH2 = multGXH([], X, IndH2);
%HYH = multGXH([], H1TY, IndH2);

% second-order part
a = multTr(KQ .* GYG .^ 2);

% first-order part
tmp1 = multTr(KQ .* GXG .* GYG);
tmp2 = multTr(GKGP .* Y);
b = 2 * tmp1 + tmp2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = optStep(a, b)
% Compute the optimal step size

t = -b / (a + eps) / 2;
if t < eps
    if a > 0
        t = 1;
    else
        t = 0;
    end
else
    if a > 0
        t = 1;
    else
        if t > 1
            t = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KP, G1, G2, H1, H2, Ct] = makeEq(KP, G1, G2, H1, H2, Ct)
% Introduce additional nodes to make the graph to be of the same size.

% dimension
[n1, m1] = size(G1);
[n2, m2] = size(G2);

if n1 < n2
    KP = [KP; zeros(n2 - n1, n2)];
    G1 = [G1; zeros(n2 - n1, m1)];
    H1 = [G1, eye(n2)];
    Ct = [Ct; ones(n2 - n1, n2)];    
else
    KP = [KP, zeros(n1, n1 - n2)];
    G2 = [G2; zeros(n1 - n2, m2)];
    H2 = [G2, eye(n1)];
    Ct = [Ct, ones(n1, n1 - n2)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, nIt, objs, objRs, X, alps)
% Debug

% score
set(gcf, 'CurrentAxes', Ax{1}); cla;
hold on;
plot(alps(1 : nIt), objs(1 : nIt), 'ro-');
plot(alps(1 : nIt), objRs(1 : nIt), 'bs-');
hT = title('objective');
set(hT, 'fontSize', 20);
set(gca, 'xlim', [0, 1]);

% x
if nIt == 1
    ha.hX = shM(X, 'ax', Ax{2});
    hT = title('correspondence');
    set(hT, 'fontSize', 20);
else
    shMUpd(ha.hX, X);
end

drawnow;
pause(1);
