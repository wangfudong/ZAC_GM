function demoHouse()
% A demo comparison of different graph matching methods on the CMU House dataset.
%
% Remark
%   The edge is directed and the edge feature is asymmetric. 
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-08-2013
warning off all;
clear variables;
clc;
close all;

addPath;
prSet(1);

dataPath = './data/cmum/graphs';
savePath = './save/cmum';

%% algorithm parameter
%[pars, algs] = gmPar(2);

for n1 = 20
    for gaps = 10:10:100
        for fid = 1:110 - gaps
            processSequence(dataPath, savePath, n1, gaps, fid);
        end
    end
end


end

function processSequence(dataPath, savePath, n1, gap, fid)


    gmat = sprintf('%s/n%d_g%d_f%d.mat', dataPath, n1, gap, fid);
    smat_ipfps = sprintf('%s/n%d_g%d_f%d_ipfps.mat', savePath, n1, gap, fid);
    smat_rrwm = sprintf('%s/n%d_g%d_f%d_rrwm.mat', savePath, n1, gap, fid);
    smat_psm = sprintf('%s/n%d_g%d_f%d_psm.mat', savePath, n1, gap, fid);
    smat_gnccp = sprintf('%s/n%d_g%d_f%d_gnccp.mat', savePath, n1, gap, fid);
    smat_bpfG = sprintf('%s/n%d_g%d_f%d_bpfG.mat', savePath, n1, gap, fid);
    smat_fgm = sprintf('%s/n%d_g%d_f%d_fgm.mat', savePath, n1, gap, fid);
    smat_gagm = sprintf('%s/n%d_g%d_f%d_gagm.mat', savePath, n1, gap, fid);
    

    n2 = 30;
    if ~exist(gmat, 'file')
        %% src parameter
        tag = 'house';
        pFs = [fid fid+gap]; % frame index
        nIn = [n1 n2];
        
        %% src
        wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
        asgT = wsSrc.asgT;

        %% feature
        parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
        parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
        wsFeat = cmumAsgFeat(wsSrc, parG, parF, 'svL', 1);
        [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

        %% affinity
        parKnl = st('alg', 'cmum'); % type of affinity: only edge distance
        [KP, KQ] = conKnlGphPQU(gphs, parKnl);
        K = conKnlGphKU(KP, KQ, gphs);
        Ct = ones(size(KP));
        [ind, ind_m] = find(Ct);

        group1 = zeros(size(ind, 1), n1);
        group2 = zeros(size(ind, 1), n2);
        for i = 1:size(ind, 1)
            group1(i, ind(i)) = 1;
            group2(i, ind_m(i)) = 1;
        end
        group1 = logical(group1);
        group2 = logical(group2);
        
        save(gmat, 'gphs', 'asgT', 'K', 'KP', 'KQ', 'group1', 'group2');
    else
        load(gmat);
        Ct = ones(size(KP));
    end

%    K = full(K);


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
    fprintf('\r\ngap = %d, t = %d, GAGM finish!', gap, fid);

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
    fprintf('\r\ngap = %d, t = %d, IPFP-S finish!', gap, fid);

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
    fprintf('\r\ngap = %d, t = %d, RRWM finish!', gap, fid);

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
    fprintf('\r\ngap = %d, t = %d, PSM finish!', gap, fid);

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
    fprintf('\r\ngap = %d, t = %d, GNCCP finish!', gap, fid);

    
    %% BPF_G
    t0 = cputime;
    parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
    X = BPF_G(K, group1, group2, parGnccp);
    X = greedyMapping(X, group1, group2);
    X = reshape(X, n1, n2);
    asgBpfG.tm = cputime - t0;
    asgBpfG.alg = 'BPF_AG';
    asgBpfG.X = X;
    asgBpfG.acc = matchAsg(X, asgT);
    asgBpfG.obj = X(:)' * K * X(:);
    fprintf('\r\ngap = %d, t = %d, BPF_G finish!', gap, fid);
    save(smat_bpfG, 'asgBpfG');


    %% FGM-U
    t0 = cputime;
    parFgmS = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n', 'ip', 'n');
    asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, parFgmS);
    asgFgmU.tm = cputime - t0;
    save(smat_fgm, 'asgFgmU');
    fprintf('\r\ngap = %d, t = %d, FgmU finish!', gap, fid);
end

