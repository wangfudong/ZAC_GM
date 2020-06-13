function demoPascal()
warning off all;
clear variables;
addPath;
prSet(1);

dataPath1 = './data/Pascal/Motorbikes';
dataPath2 = './data/Pascal/Cars';
graphPath = './data/Pascal/Graphs';
savePath = './save/Pascal';

outliers = 0:2:20;
for cImg =1:50
    ProcessImage(dataPath1, dataPath2, graphPath, savePath, cImg, outliers);  
end
end


function ProcessImage(dataPath1, dataPath2, graphPath, savePath, cImg, outliers)
if cImg <= 20
    fmat = sprintf('%s/pair_%d.mat', dataPath1, cImg);
    gmat = sprintf('%s/Motor%02d.mat', graphPath, cImg);
    smat_ipfps = sprintf('%s/Motor%02d_ipfps.mat', savePath, cImg);
    smat_psm = sprintf('%s/Motor%02d_psm.mat', savePath, cImg);
    smat_rrwm = sprintf('%s/Motor%02d_rrwm.mat', savePath, cImg);
    smat_gnccp = sprintf('%s/Motor%02d_gnccp.mat', savePath, cImg);
    smat_bpfG = sprintf('%s/Motor%02d_bpfG.mat', savePath, cImg);
    smat_fgm = sprintf('%s/Motor%02d_fgm.mat', savePath, cImg);
    smat_gagm = sprintf('%s/Motor%02d_gagm.mat', savePath, cImg);
else
    fmat = sprintf('%s/pair_%d.mat', dataPath2, cImg - 20);
    gmat = sprintf('%s/Car%02d.mat', graphPath, cImg - 20);
    smat_ipfps = sprintf('%s/Car%02d_ipfps.mat', savePath, cImg - 20);
    smat_psm = sprintf('%s/Car%02d_psm.mat', savePath, cImg - 20);
    smat_rrwm = sprintf('%s/Car%02d_rrwm.mat', savePath, cImg - 20);
    smat_gnccp = sprintf('%s/Car%02d_gnccp.mat', savePath, cImg - 20);
    smat_bpfG = sprintf('%s/Car%02d_bpfG.mat', savePath, cImg - 20);%
    smat_fgm = sprintf('%s/Car%02d_fgm.mat', savePath, cImg - 20);
    smat_gagm = sprintf('%s/Car%02d_gagm.mat', savePath, cImg - 20);
end
    

if exist(gmat, 'file')
   load(gmat);
else
    GRAPH = cell(size(outliers));
end

asgIpfpS = cell(size(outliers));
asgRrwm = cell(size(outliers));
asgPsm = cell(size(outliers));
asgGnccp = cell(size(outliers));
asgBpfG = cell(size(outliers));
asgFgm = cell(size(outliers));
asgGagm = cell(size(outliers));


   for i = 1:size(outliers, 2)
        if ~exist(gmat, 'file')
            [K, group1, group2, KP, KQ, gphs, X_GT] = conPascalGphs(fmat, outliers(i));
            GRAPH{i}.K = K;
            GRAPH{i}.KP = KP;
            GRAPH{i}.KQ = KQ;
            GRAPH{i}.group1 = group1;
            GRAPH{i}.group2 = group2;
            GRAPH{i}.gphs = gphs;
            GRAPH{i}.X_GT = X_GT;
        else
            K = GRAPH{i}.K;
            KP = GRAPH{i}.KP;
            KQ = GRAPH{i}.KQ;
            group1 = GRAPH{i}.group1;
            group2 = GRAPH{i}.group2;
            gphs = GRAPH{i}.gphs;
            X_GT = GRAPH{i}.X_GT;
        end
        asgT.X = X_GT;
     
        K = sparse(K);

        Ct = ones(size(KP));
        [n1, n2] = size(KP);

        
        %% GAGM
        t0 = cputime;
        E12 = ones(size(group1, 2), size(group2,2));
        [L1, L2] = find(E12);
        L12 = [L1, L2];
        X = GAM(K, L12, group1, group2);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgGagm{i}.tm = cputime - t0;
        asgGagm{i}.alg = 'GAGM';
        asgGagm{i}.X = X;
        asgGagm{i}.acc = matchAsg(X, asgT);
        asgGagm{i}.obj =  X(:)' * K * X(:);
         fprintf('\r\ncImg = %d, outliers = %d, GAGM finish!', cImg, outliers(i));
        
        %% IPFP-S
        t0 = cputime;
        X = IPFP_S(K, group1, group2);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgIpfpS{i}.tm = cputime - t0;
        asgIpfpS{i}.alg = 'IpfpS';
        asgIpfpS{i}.X = X;
        asgIpfpS{i}.acc = matchAsg(X, asgT);
        asgIpfpS{i}.obj = X(:)' * K * X(:);
        fprintf('\r\ncImg = %d, outliers = %d, IPFP_S finish!', cImg, outliers(i));
        
        %% RRWM
        t0 = cputime;
        X = RRWM(K, group1, group2);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgRrwm{i}.tm = cputime - t0;
        asgRrwm{i}.alg = 'rrwm';
        asgRrwm{i}.X = X;
        asgRrwm{i}.acc = matchAsg(X, asgT);
        asgRrwm{i}.obj = X(:)' * K * X(:);
        fprintf('\r\ncImg = %d, outliers = %d, RRWM finish!', cImg, outliers(i));

        %% PSM
        t0 = cputime;
        X = PSM(K, group1, group2);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgPsm{i}.tm = cputime - t0;
        asgPsm{i}.alg = 'psm';
        asgPsm{i}.X = X;
        asgPsm{i}.acc = matchAsg(X, asgT);
        asgPsm{i}.obj = X(:)' * K * X(:);
        fprintf('\r\ncImg = %d, outliers = %d, PSM finish!', cImg, outliers(i));
        
        %% GNCCP
        t0 = cputime;
        parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
        X = GNCCP_K(K, group1, group2, parGnccp);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgGnccp{i}.tm = cputime - t0;
        asgGnccp{i}.alg = 'gnccp';
        asgGnccp{i}.X = X;
        asgGnccp{i}.acc = matchAsg(X, asgT);
        asgGnccp{i}.obj = X(:)' * K * X(:);
        fprintf('\r\ncImg = %d, outliers = %d, GNCCP finish!', cImg, outliers(i));
    
        
        %% BPF_G
        t0 = cputime;
        parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
        X = BPF_G(K, group1, group2, parGnccp);
        X = greedyMapping(X, group1, group2);
        X = reshape(X, n1, n2);
        asgBpfG{i}.tm = cputime - t0;
        asgBpfG{i}.alg = 'BPF_AG';
        asgBpfG{i}.X = X;
        asgBpfG{i}.acc = matchAsg(X, asgT);
        asgBpfG{i}.obj = X(:)' * K * X(:);
        fprintf('\r\ncImg = %d, outliers = %d, BPF_G finish!\r\n', cImg, outliers(i));

        %% FGM-U
        t0 = cputime;
        parFgmS = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n', 'ip', 'n');
        asgFgm{i} = fgmU(KP, KQ, Ct, gphs, asgT, parFgmS);
        asgFgm{i}.tm = cputime - t0;
        fprintf('\r\ncImg = %d, outliers = %d, FGM finish!', cImg, outliers(i));

   end
 
if ~exist(gmat, 'file')
   save(gmat, 'GRAPH'); 
end
save(smat_ipfps, 'asgIpfpS');
save(smat_rrwm, 'asgRrwm');
save(smat_psm, 'asgPsm');
save(smat_gnccp, 'asgGnccp');
save(smat_bpfG, 'asgBpfG');
save(smat_fgm, 'asgFgm');
save(smat_gagm, 'asgGagm');
end

