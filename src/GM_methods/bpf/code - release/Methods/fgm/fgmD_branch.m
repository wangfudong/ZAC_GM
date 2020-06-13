function asg = fgmD_branch(KP, KQ, Ct, gphs, asgT, par)
% Factorized graph matching.
%
% Remark
%   The edge is directed and the edge feature is asymmetric.
%
% Reference
%   F. Zhou and F. De la Torre, "Deformable Graph Matching", in CVPR, 2013.
%
% Input
%   KP       -  node affinity matrix, n1 x n2
%   KQ       -  edge affinity matrix, m1 x m2
%   Ct       -  correspondence constraint, n1 x n2
%                 Ct_ij = 1: i and j can be matched
%                 Ct_ij = 0: i and j cannot be matched
%   gphs     -  graphs, 1 x 2 (cell)
%     G      -  node-edge adjacency (for starting point), n x mi
%     H      -  node-edge adjacency (for ending point), n x mi
%   par      -  parameter
%     nAlp   -  #alpha or #steps to do the path-following, {100}
%     nItMa  -  #maximum iteration steps for each scale of alpha, {100}
%     nHst   -  #history nodes for modifed FW algorithm, {10}
%     ip     -  flag of using IPFP to improve the algorithm, {'y'} | 'n'
%     deb    -  flag of debugging, 'y' | {'n'}
%
% Output
%   asg      -  assignment
%     alg    -  'fgmD'
%     X      -  correspondence matrix, n1 x n2
%     acc    -  accuracy
%     obj    -  objective
%
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-08-2013

% function parameter
nAlp = ps(par, 'nAlp', 101);
nItMa = ps(par, 'nItMa', 100);
nHst = ps(par, 'nHst', 10);
isIp = psY(par, 'ip', 'n');
isDeb = psY(par, 'deb', 'n');
prIn('fgmD', 'nAlp %d, nItMa %d, nHst %d, isIp %d', ...
     nAlp, nItMa, nHst, isIp);

% weight
alps = linspace(0, 1, nAlp);

% store variables
XQ1 = ps(par, 'XQ1', []);
XQ2 = ps(par, 'XQ2', []);
pathDStore(KP, KQ, Ct, XQ1, XQ2, gphs, 'path', []);


% original path-following
X0 = gmIniUnif(Ct, st('nor', 'doub'));
initAlp = 1;
[path, X, ~, obj, nIts, Xs, objs, objGms, objCons, objVexs, objCavs, useIps, objInss, objIn2ss] = ...
    pathDIter(X0, initAlp, alps, nItMa, nHst, isIp, isDeb, 1, [], 0);
optX = X;
maxObj = obj;
count  = 1;

g_pre = 0;
g_err = zeros(1, nAlp);
for i = 2 : size(path, 1)
    g = reshape(path(i,:,:) - path(i-1, :, :), size(path, 2) * size(path, 3), 1);
    g_err(i) = norm(g - g_pre);
    g_pre = g;
    
    if i <= nHist
        continue;
    end
    
    if ~isSingularPoint(g_err(1:i), nHist)
        continue;
    end
    
    fprintf('\n  === iAlp = %d, Singular point detected, branch a new path!!! ===', i - 1);
    
     steps = 1;
     alp_arc = alps(i-1+steps);
     X_arc = arclength(path, i-1, steps);
     X_arc = bistocNormalize_slack(X_arc, 1e-3);
     
        
    [~, X, ~, obj, nIts, Xs, objs, objGms, objCons, objVexs, objCavs, useIps, objInss, objIn2ss] = ...
        pathDIter(X_arc, alp_arc, alps, nItMa, nHst, isIp, isDeb, 1, [], 0);
    if obj > maxObj
        maxObj = obj;
        optX = X;
    end
    
    count = count + 1;
end

fprintf('\n ** Totaolly %d branches are searched! ** \n', count);

X = optX;
obj = maxObj;

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% store
asg.alg = 'fgmD';
asg.X = X;
asg.obj = full(obj);
asg.acc = acc;

% for debug
asg.nIts = nIts;
asg.Xs = Xs;
asg.objs = objs;
asg.objGms = objGms;
asg.objVexs = objVexs;
asg.objCavs = objCavs;
asg.objCons = objCons;
asg.objInss = objInss;
asg.objIn2ss = objIn2ss;
asg.useIps = useIps;

prOut;

function X_arc = arclength(path, ind, steps, nHist)
x_len = size(path, 2) * size(path, 3);
x_store = zeros(x_len, 1);
for i = ind - nHist : ind
    x = reshape(path(i, :, :), x_len, 1);
    x_store = [x_store, x];
end
x_gra = EstimateGradient(x_store);
x_arc = x + x_gra * steps;
X_arc = reshape(x_arc, size(path, 2), size(path, 3));


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

