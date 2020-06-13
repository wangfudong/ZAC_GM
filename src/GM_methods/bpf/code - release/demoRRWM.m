function demoRRWM()
% A demo comparison of different graph matching methods on the dataset provided by RRWM .
%
% History
%   create  -  Tao WANG (twang@bjtu.edu.cn), 01-02-2015
clear variables;
%prSet(1);

dataPath = './data/RRWM';  % Path for 'mat' files
savePath = './save/RRWM/';


for imgIndex = 1:30
    processImg(dataPath, savePath, imgIndex);
end

end

function processImg(dataPath, savePath, imgIndex)
     fprintf('\n\n******** image pair %02d *****\n', imgIndex);   
     
    %% load data
    datafile = sprintf('%s/matchData/fi_%02da+%02db.mat', dataPath, imgIndex, imgIndex);
    load(datafile);

    %% Options & parameters for experiment
    affinity_max = 50;               % maximum value of affinity 
    
    %% affinity matrix
    % Make affinity matrix
    cdata.affinityMatrix = max(affinity_max - cdata.distanceMatrix,0); % dissimilarity -> similarity conversion
    %cdata.affinityMatrix = exp(-cdata.distanceMatrix/25); % dissimilarity -> similarity conversion
    cdata.affinityMatrix(1:(length(cdata.affinityMatrix)+1):end) = 0; % diagonal zeros

    K = cdata.affinityMatrix;
    group1 = full(cdata.group1);
    group2 = full(cdata.group2);
    
    savefile = sprintf('%s%02d_MatchX.mat', savePath, imgIndex); 
    
    
   
    t0 = cputime;
    MatchX.IPFP = IPFP_S(K, group1, group2);
    MatchX.IPFP = greedyMapping(MatchX.IPFP, group1, group2);
    CT.IPFP = cputime - t0;
    fprintf('\ncImg = %d, IPFPS finished!', imgIndex);
    
    t0 = cputime;
    MatchX.RRWM = RRWM(K, group1, group2);
    MatchX.RRWM = greedyMapping(MatchX.RRWM, group1, group2);
    CT.RRWM = cputime - t0;
    fprintf('\ncImg = %d, RRWM finished!', imgIndex);
    
    t0 = cputime;
    MatchX.PSM = PSM(K, group1, group2);
    MatchX.PSM = greedyMapping(MatchX.PSM, group1, group2);
    CT.PSM = cputime - t0;
    fprintf('\ncImg = %d, PSM finished!', imgIndex);

    t0 = cputime;
    parGnccp = st('nItMa', 100, 'deta', 0.01, 'nHist', 5, 'rho', 2, 'theta', 0.01);
    MatchX.GNCCP = GNCCP_K(K, group1, group2, parGnccp);
    MatchX.GNCCP = greedyMapping(MatchX.GNCCP, group1, group2);
    CT.GNCCP = cputime-t0;
    fprintf('\ncImg = %d, GNCCP finished!', imgIndex);

    
    t0 = cputime;
    parGnccp = st('nItMa', 100, 'deta', 0.01, 'nHist', 5, 'rho', 2, 'theta', 0.05);
    MatchX.BPF_G = BPF_G(K, group1, group2, parGnccp);
    MatchX.BPF_G = greedyMapping(MatchX.BPF_G, group1, group2);
    CT.BPF_G = cputime - t0;
    fprintf('\ncImg = %d, BPF_G finished!', imgIndex);

    save(savefile, 'MatchX', 'CT');        
end
