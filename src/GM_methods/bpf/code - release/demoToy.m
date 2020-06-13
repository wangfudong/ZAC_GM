% A demo comparison of different graph matching methods on the synthetic dataset.
%
% Remark
%   The edge is directed and the edge feature is asymmetric.
%

function demoToy
warning off all;
clear variables;
addPath;
prSet(1);


    setting = 'outliers';
    savePath = ['./save/Toy/' setting];
    egDen = .5; % edge density
    egDef = 0; % edge deformation
    for outliers = 0:10
        for gid = 1:100
            processSyntheticGph(savePath, setting, outliers, egDen, egDef, gid);
        end
    end
    
    setting = 'noise';
    savePath = ['./save/Toy/' setting];
    egDen = .5;
    outliers = 0;
    for egDef = 0:0.02:0.2
        for gid = 1:100
            processSyntheticGph(savePath, setting, outliers, egDen, egDef, gid);            
        end
    end

   setting = 'density';
   savePath = ['./save/Toy/' setting];
    outliers = 0;
    egDef = 0.1;
    for egDen = 0.3:0.1:1
        for gid = 1:100
            processSyntheticGph(savePath, setting, outliers, egDen, egDef, gid);            
        end
    end

end

function processSyntheticGph(savePath, setting, outliers, egDen, egDef, gid)

if strcmp(setting, 'outliers')
    var = outliers;
elseif strcmp(setting, 'noise')
    var = round(egDef * 100);
elseif strcmp(setting, 'density')
    var = round(egDen * 10);
end

    smat_ipfps  = sprintf('%s/%s%02d_%d_ipfps.mat', savePath, setting, var, gid);
    smat_rrwm   = sprintf('%s/%s%02d_%d_rrwm.mat', savePath, setting, var, gid);
    smat_psm    = sprintf('%s/%s%02d_%d_psm.mat', savePath, setting, var, gid);
    smat_gnccp  = sprintf('%s/%s%02d_%d_gnccp.mat', savePath, setting, var, gid);
    smat_BpfG   = sprintf('%s/%s%02d_%d_BpfG.mat', savePath,setting, var, gid);
    smat_fgm   = sprintf('%s/%s%02d_%d_fgm.mat', savePath,setting, var, gid);
    smat_gagm   = sprintf('%s/%s%02d_%d_gagm.mat', savePath,setting, var, gid);

    %% src parameter
    tag = 1; 
    nIn = 20; % #inliers
    nOuts = [outliers outliers]; % #outliers
    parKnl = st('alg', 'toy'); % type of affinity: synthetic data

    %% src
    wsSrc = toyAsgSrcD(tag, nIn, nOuts, egDen, egDef);
    [gphs, asgT] = stFld(wsSrc, 'gphs', 'asgT');

    %% affinity
    [KP, KQ] = conKnlGphPQD(gphs, parKnl); % node and edge affinity
    K = conKnlGphKD(KP, KQ, gphs); % global affinity
    Ct = ones(size(KP)); % mapping constraint (default to a constant matrix of one)

%     % directed graph -> undirected graph (for fgmU)
%     gphUs = gphD2Us(gphs);
%     [~, KQU] = knlGphKD2U(KP, KQ, gphUs);

    % group info
    [ind, ind_m] = find(Ct);
    n1 = nIn + nOuts(1);
    n2 = nIn + nOuts(2);
    group1 = zeros(size(ind, 1), n1);
    group2 = zeros(size(ind, 1), n2);
    for i = 1:size(ind, 1)
        group1(i, ind(i)) = 1;
        group2(i, ind_m(i)) = 1;
    end
    group1 = logical(group1);
    group2 = logical(group2);

    %% GAGM
    t0 = cputime;
    E12 = ones(size(group1, 2), size(group2,2));
    [L1, L2] = find(E12);
    L12 = [L1, L2];
    X = GAM(K, L12, group1, group2);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgGagm.tm = cputime - t0;
    asgGagm.alg = 'GAGM';
    asgGagm.X = X;
    asgGagm.acc = matchAsg(X, asgT);
    asgGagm.obj = X(:)' * K * X(:);
    save(smat_gagm, 'asgGagm');
    fprintf('\r\n %s = %d, t = %d, GAGM finish!', setting, var, gid);
    
    %% IPFP-S
    t0 = cputime;
    X = IPFP_S(K, group1, group2);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgIpfpS.tm = cputime - t0;
    asgIpfpS.alg = 'IpfpS';
    asgIpfpS.X = X;
    asgIpfpS.acc = matchAsg(X, asgT);
    asgIpfpS.obj = X(:)' * K * X(:);
    save(smat_ipfps, 'asgIpfpS');
    fprintf('\r\n %s = %d, t = %d, IPFP-S finish!', setting, var, gid);

    %% RRWM
    t0 = cputime;
    X = RRWM(K, group1, group2);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgRrwm.tm = cputime - t0;
    asgRrwm.alg = 'rrwm';
    asgRrwm.X = X;
    asgRrwm.acc = matchAsg(X, asgT);
    asgRrwm.obj = X(:)' * K * X(:);
    save(smat_rrwm, 'asgRrwm');
    fprintf('\r\n %s = %d, t = %d, RRWM finish!', setting, var, gid);

    %% PSM
    t0 = cputime;
    X = PSM(K, group1, group2);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgPsm.tm = cputime - t0;
    asgPsm.alg = 'psm';
    asgPsm.X = X;
    asgPsm.acc = matchAsg(X, asgT);
    asgPsm.obj = X(:)' * K * X(:);
    save(smat_psm, 'asgPsm');
    fprintf('\r\n %s = %d, t = %d, PSM finish!', setting, var, gid);

  
    %% GNCCP
    t0 = cputime;
    parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
    X = GNCCP_K(K, group1, group2, parGnccp);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgGnccp.tm = cputime - t0;
    asgGnccp.alg = 'gnccp';
    asgGnccp.X = X;
    asgGnccp.acc = matchAsg(X, asgT);
    asgGnccp.obj = X(:)' * K * X(:);
    save(smat_gnccp, 'asgGnccp');
    fprintf('\r\n %s = %d, t = %d, GNCCP finish!', setting, var, gid);

    %% BPF_G
    t0 = cputime;
    parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
    X = BPF_G(K, group1, group2, parGnccp);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgBpfG.tm = cputime - t0;
    asgBpfG.alg = 'BPF_G';
    asgBpfG.X = X;
    asgBpfG.acc = matchAsg(X, asgT);
    asgBpfG.obj = X(:)' * K * X(:);
    save(smat_BpfG, 'asgBpfG');
    fprintf('\r\n %s = %d, t = %d, BPF_G finish!', setting, var, gid);

   
    %% FGM-D
    t0 = cputime;
    parFgmS = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n', 'ip', 'n');
    asgFgmD = fgmD(KP, KQ, Ct, gphs, asgT, parFgmS);
    asgFgmD.tm = cputime - t0;
    save(smat_fgm, 'asgFgmD');
    fprintf('\r\n %s = %d, t = %d, FGM finish!', setting, var, gid);
end