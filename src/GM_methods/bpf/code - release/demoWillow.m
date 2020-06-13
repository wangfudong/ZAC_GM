function demoWillow()
warning off all;
clear variables;
close all;
addPath;
prSet(1);

imgPath  = './data/Willow';
dataPath = './data/Willow/graphs';
savePath = './save/Willow';
imgClass = {'Car', 'Duck', 'Face', 'Motorbike', 'Winebottle'};

for cls = 1 : size(imgClass, 2)
    for cImg = 1:100
        ProcessClass(imgPath, dataPath, savePath, imgClass{cls}, cImg);
    end
end

end

function ProcessClass(imgPath, dataPath, savePath, clsName, cImg)
searchDir = [imgPath '/' clsName '/*.png'];    
lstImgs = dir(searchDir);    
nImages = size(lstImgs, 1);

outliers = 0:10;

asgIpfpS = cell(size(outliers));
asgRrwm = cell(size(outliers));
asgPsm = cell(size(outliers));
asgGnccp = cell(size(outliers));
asgBpfG = cell(size(outliers));
asgFgmD = cell(size(outliers));
asgGagm = cell(size(outliers));

    dataFile = sprintf('%s/%s_%04d.mat', dataPath, clsName, cImg);
    
    smat_ipfps  = sprintf('%s/%s_%04d_ipfpS.mat', savePath, clsName, cImg);
    smat_psm    = sprintf('%s/%s_%04d_psm.mat', savePath, clsName, cImg);
    smat_rrwm   = sprintf('%s/%s_%04d_rrwm.mat', savePath, clsName, cImg);
    smat_gnccp  = sprintf('%s/%s_%04d_gnccp.mat', savePath, clsName, cImg);
    smat_BpfG = sprintf('%s/%s_%04d_BpfG.mat', savePath, clsName, cImg);
    smat_fgm = sprintf('%s/%s_%04d_fgm.mat', savePath, clsName, cImg);
    smat_gagm = sprintf('%s/%s_%04d_gagm.mat', savePath, clsName, cImg);
    
    fprintf('\n cImg = %d, start processing!', cImg);
    
    if exist(dataFile, 'file')
        fprintf('\n cImg = %d, data File exist!!!', cImg);
        load(dataFile);
    else
        index = randperm(nImages, 2);
        
        img1  = [imgPath '/' clsName '/' lstImgs(index(1)).name];
        anno1 = [imgPath '/' clsName '/' lstImgs(index(1)).name(1:end-4) '.mat'];
        img2  = [imgPath '/' clsName '/' lstImgs(index(2)).name];
        anno2 = [imgPath '/' clsName '/' lstImgs(index(2)).name(1:end-4) '.mat'];

        clear GRAPH;
        GRAPH = cell(size(outliers));
        for i = 1:size(outliers, 2)
            GRAPH{i}.cdata = conWillowGphs(img1, anno1, img2, anno2, outliers(i));
            GRAPH{i}.nOutliers = outliers(i);
        end
        
        save(dataFile, 'GRAPH', 'index');
    end
    
    for i = 1:size(outliers, 2)
        cdata = GRAPH{i}.cdata;
        n1 = size(cdata.group1, 2);
        n2 = size(cdata.group2, 2);
        asgT.X = cdata.X_GT;
        Ct = ones(n1, n2);
        
        % set node affinity to zero;
     %   cdata.K(1:(size(cdata.K,1)+1):end) = 0;
     %   avgv = sum(cdata.K(:)) / sum(cdata.K(:)>0);
     %   cdata.K = max(0, cdata.K - avgv);

        cdata.K = cdata.KD;
        cdata.KQ = cdata.KQD;
        cdata.gphs = cdata.gphDs;

        %% GAGM
        t0 = cputime;
        E12 = ones(size(cdata.group1, 2), size(cdata.group2,2));
        [L1, L2] = find(E12);
        L12 = [L1, L2];
        x = GAM(cdata.K, L12, cdata.group1, cdata.group2);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgGagm{i}.tm = cputime - t0;
        asgGagm{i}.alg = 'GAGM';
        asgGagm{i}.X = X;
        asgGagm{i}.acc = matchAsg(X, asgT);
        asgGagm{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, GAGM finish!', cImg,index(1), index(2), outliers(i));
        
        %% IPFP-S
        t0 = cputime;
        x = IPFP_S(cdata.K, cdata.group1, cdata.group2);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgIpfpS{i}.tm = cputime - t0;
        asgIpfpS{i}.alg = 'IpfpS';
        asgIpfpS{i}.X = X;
        asgIpfpS{i}.acc = matchAsg(X, asgT);
        asgIpfpS{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, IPFP_S finish!', cImg,index(1), index(2), outliers(i));
        
        %% RRWM
        t0 = cputime;
        x = RRWM(cdata.K, cdata.group1, cdata.group2);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgRrwm{i}.tm = cputime - t0;
        asgRrwm{i}.alg = 'rrwm';
        asgRrwm{i}.X = X;
        asgRrwm{i}.acc = matchAsg(X, asgT);
        asgRrwm{i}.obj = x' * cdata.K *x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, RRWM finish!', cImg, index(1), index(2),outliers(i));

        %% PSM
        t0 = cputime;
        x = PSM(cdata.K, cdata.group1, cdata.group2);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgPsm{i}.tm = cputime - t0;
        asgPsm{i}.alg = 'psm';
        asgPsm{i}.X = X;
        asgPsm{i}.acc = matchAsg(X, asgT);
        asgPsm{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, PSM finish!', cImg, index(1), index(2),outliers(i));
        
        %% GNCCP
        t0 = cputime;
        parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
        x = GNCCP_K(cdata.K, cdata.group1, cdata.group2, parGnccp);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgGnccp{i}.tm = cputime - t0;
        asgGnccp{i}.alg = 'gnccp';
        asgGnccp{i}.X = X;
        asgGnccp{i}.acc = matchAsg(X, asgT);
        asgGnccp{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, GNCCP finish!', cImg,index(1), index(2), outliers(i));
        
        %% BPF_G
        t0 = cputime;
        parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
        x = BPF_G(cdata.K, cdata.group1, cdata.group2, parGnccp);
        x = greedyMapping(x, cdata.group1, cdata.group2);
        X = TransToMatrix(x, cdata.matchInfo, n1, n2);
        asgBpfG{i}.tm = cputime - t0;
        asgBpfG{i}.alg = 'BPF_G';
        asgBpfG{i}.X = X;
        asgBpfG{i}.acc = matchAsg(X, asgT);
        asgBpfG{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, BPF_G finish!', cImg,index(1), index(2), outliers(i));

        %% FGM-D
        t0 = cputime;
        parFgmS = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n', 'ip', 'n');
        asgFgmD{i} = fgmD(cdata.KP, cdata.KQD, Ct, cdata.gphDs, asgT, parFgmS);
        x = TransToVector(asgFgmD{i}.X, cdata.matchInfo, n1, n2);
        asgFgmD{i}.tm = cputime - t0;
        asgFgmD{i}.obj = x' * cdata.K * x;
        fprintf('\r\ncImg = %d(%d,%d), outliers = %d, FGM finish!', cImg, index(1), index(2),outliers(i));
           
    end
    
   save(smat_ipfps, 'asgIpfpS');
   save(smat_rrwm, 'asgRrwm');
   save(smat_psm, 'asgPsm');
   save(smat_gnccp, 'asgGnccp');
   save(smat_BpfG, 'asgBpfG');
   save(smat_fgm, 'asgFgmD');
   save(smat_gagm, 'asgGagm');
end

function X = TransToMatrix(x, matchInfo, n1, n2)
X = zeros(n1, n2);
for i = 1 : n1 * n2
    X(matchInfo(1,i), matchInfo(2,i)) = x(i);
end
end

function x = TransToVector(X, matchInfo, n1, n2)
x = zeros(n1 * n2, 1);
for i = 1 : n1 * n2
    x(i) = X(matchInfo(1,i), matchInfo(2,i));
end
end
