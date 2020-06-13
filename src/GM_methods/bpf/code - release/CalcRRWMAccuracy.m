function CalcRRWMAccuracy()
dataPath = './data/RRWM';  % Path for 'mat' files
savePath = './save/RRWM/';
fid = fopen('RRWM.txt', 'w');

%% Options & parameters for experiment
bDisplayMatching = 0;            % Display image feature matching results or not
extrapolate_thres = 15;          % Extrapolation of matches for flexible evaluation
affinity_max = 50;               % maximum value of affinity 

avgAcc.Ipfp = 0;
avgAcc.Rrwm = 0;
avgAcc.Psm = 0;
avgAcc.Gnccp = 0;
avgAcc.BpfG = 0; %zeros(size(nHist, 2), size(theta,2));


avgCT = avgAcc;

curAcc = avgAcc;
objectives = avgAcc;
curObj = objectives;

bRefineX = 1;


for imgIndex = 1:30
    
    fprintf('\n\n******** image pair %02d *****\n', imgIndex);   

    %% load data
    datafile = sprintf('%s/matchData/fi_%02da+%02db.mat', dataPath, imgIndex, imgIndex);
    xfile = sprintf('%s%02d_MatchX.mat', savePath, imgIndex); 
    load(datafile);
    load(xfile);

    %% affinity matrix
    cdata.affinityMatrix = max(affinity_max - cdata.distanceMatrix,0); % dissimilarity -> similarity conversion
    %cdata.affinityMatrix = exp(-cdata.distanceMatrix/25); % dissimilarity -> similarity conversion
    cdata.affinityMatrix(1:(length(cdata.affinityMatrix)+1):end) = 0; % diagonal zeros
    K = cdata.affinityMatrix;
    
    % Extrapolate the given ground truths for flexible evaluation
    cdata.GTbool = extrapolateGT(cdata.view, cell2mat({cdata.matchInfo.match}'), cdata.GT, extrapolate_thres)';
    cdata.extrapolate_thres = extrapolate_thres;
    
    X = MatchX.IPFP;
    if bRefineX
        X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
    end
    curObj.Ipfp = MatchX.IPFP'*K*MatchX.IPFP;
    objectives.Ipfp = objectives.Ipfp + curObj.Ipfp;
    curAcc.Ipfp = sum(cdata.GTbool .* X) / sum(cdata.GTbool);
    avgAcc.Ipfp = avgAcc.Ipfp + curAcc.Ipfp;
    avgCT.Ipfp = avgCT.Ipfp + CT.IPFP;
    if bDisplayMatching
        str = ['cImg=' num2str(imgIndex) ',IPFP,' ...
            num2str(sum(cdata.GTbool .* X)) '/' num2str(sum(cdata.GTbool))];
        figure('NumberTitle', 'off', 'Name', str);
        displayFeatureMatching(cdata, X, cdata.GTbool);
    end
    
    X = MatchX.RRWM;
    if bRefineX
        X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
    end
    curObj.Rrwm = MatchX.RRWM'*K*MatchX.RRWM;
    objectives.Rrwm = objectives.Rrwm + curObj.Rrwm;
    curAcc.Rrwm = sum(cdata.GTbool .* X) / sum(cdata.GTbool);
    avgAcc.Rrwm = avgAcc.Rrwm + curAcc.Rrwm;
    avgCT.Rrwm = avgCT.Rrwm + CT.RRWM;
    if bDisplayMatching
        str = ['cImg=' num2str(imgIndex) ',RRWM,' ...
            num2str(sum(cdata.GTbool .* X))  '/' num2str(sum(cdata.GTbool))];
        figure('NumberTitle', 'off', 'Name', str);
        displayFeatureMatching(cdata, X, cdata.GTbool);
    end

    X = MatchX.PSM;
    if bRefineX
        X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
    end
    curObj.Psm = MatchX.PSM'*K*MatchX.PSM;
    objectives.Psm = objectives.Psm + curObj.Psm;
    curAcc.Psm = sum(cdata.GTbool .* X) / sum(cdata.GTbool);
    avgAcc.Psm = avgAcc.Psm + curAcc.Psm;
    avgCT.Psm = avgCT.Psm + CT.PSM;
    if bDisplayMatching
        str = ['cImg=' num2str(imgIndex) ',PSM,' ...
            num2str(sum(cdata.GTbool .* X)) '/' num2str(sum(cdata.GTbool))];
        figure('NumberTitle', 'off', 'Name', str);
        displayFeatureMatching(cdata, X, cdata.GTbool);
    end

            
    X = MatchX.GNCCP;
    if bRefineX
        X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
    end
    curObj.Gnccp = MatchX.GNCCP'*K*MatchX.GNCCP;
    objectives.Gnccp = objectives.Gnccp + curObj.Gnccp;
    curAcc.Gnccp = sum(cdata.GTbool .* X) / sum(cdata.GTbool);
    avgAcc.Gnccp = avgAcc.Gnccp + curAcc.Gnccp;
    avgCT.Gnccp = avgCT.Gnccp + CT.GNCCP;
    if bDisplayMatching
        str = ['cImg=' num2str(imgIndex) ',GNCCP,' ...
            num2str(sum(cdata.GTbool .* X)) '/' num2str(sum(cdata.GTbool))];
        figure('NumberTitle', 'off', 'Name', str);
        displayFeatureMatching(cdata, X, cdata.GTbool);
    end
    
    X = MatchX.BPF_G;
    if bRefineX
        X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
    end
    curObj.BpfG = MatchX.BPF_G'*K*MatchX.BPF_G;
    objectives.BpfG = objectives.BpfG +  curObj.BpfG;
    curAcc.BpfG = sum(cdata.GTbool .* X) / sum(cdata.GTbool);
    avgAcc.BpfG = avgAcc.BpfG + curAcc.BpfG;
    avgCT.BpfG = avgCT.BpfG + CT.BPF_G;
    if bDisplayMatching
        str = ['cImg=' num2str(imgIndex) ',BPF_G,' ...
            num2str(sum(cdata.GTbool .* X)) '/' num2str(sum(cdata.GTbool))];
        figure('NumberTitle', 'off', 'Name', str);
        displayFeatureMatching(cdata, X, cdata.GTbool);
    end    

        
    %% show current accuracy
    fprintf(fid, ' \r\n============== algorithm Accuracy, image pair %d ======== ', imgIndex);
    fprintf(fid, '\r\nIPFP    : acc: %.2f,  objectives: %.2f', curAcc.Ipfp,     curObj.Ipfp);
    fprintf(fid, '\r\nRRWM    : acc: %.2f,  objectives: %.2f', curAcc.Rrwm,     curObj.Rrwm);
    fprintf(fid, '\r\nPSM     : acc: %.2f,  objectives: %.2f', curAcc.Psm,      curObj.Psm);
    fprintf(fid, '\r\nGNCCP   : acc: %.2f,  objectives: %.2f', curAcc.Gnccp,    curObj.Gnccp);
    fprintf(fid, '\r\nBPF_G: acc: %.2f,  objectives: %.2f',     curAcc.BpfG,    curObj.BpfG);

    fprintf(' \r\n============== algorithm Accuracy, image pair %d ======== ', imgIndex);
    fprintf('\r\nIPFP    : acc: %.2f,  objectives: %.2f', curAcc.Ipfp,    curObj.Ipfp);
    fprintf('\r\nRRWM    : acc: %.2f,  objectives: %.2f', curAcc.Rrwm,     curObj.Rrwm);
    fprintf('\r\nPSM     : acc: %.2f,  objectives: %.2f', curAcc.Psm,      curObj.Psm);
    fprintf('\r\nGNCCP   : acc: %.2f,  objectives: %.2f', curAcc.Gnccp,    curObj.Gnccp);
    fprintf('\r\nBPF_G   : acc: %.2f,  objectives: %.2f', curAcc.BpfG,      curObj.BpfG);
    
    
end

t = imgIndex;

max_obj = max([objectives.Ipfp, objectives.Rrwm, objectives.Psm, objectives.Gnccp, objectives.BpfG]);

%% show average accuracy
fprintf('\r\n============== Average Accuracy =================== ');
fprintf('\r\nIPFP    : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Ipfp / t,  objectives.Ipfp / max_obj, avgCT.Ipfp / t);
fprintf('\r\nRRWM    : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Rrwm / t,   objectives.Rrwm / max_obj,  avgCT.Rrwm / t);
fprintf('\r\nPSM     : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Psm / t,    objectives.Psm / max_obj,   avgCT.Psm / t);
fprintf('\r\nGNCCP   : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Gnccp / t,  objectives.Gnccp / max_obj, avgCT.Gnccp / t);
fprintf('\r\nBPF_G   : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.BpfG / t,  objectives.BpfG / max_obj, avgCT.BpfG / t);


%% show average accuracy
fprintf(fid, '\r\n============== Average Accuracy =================== ');
fprintf(fid,'\r\nIPF     : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Ipfp / t,  objectives.Ipfp / max_obj, avgCT.Ipfp / t);
fprintf(fid,'\r\nRRWM    : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Rrwm / t,   objectives.Rrwm / max_obj,  avgCT.Rrwm / t);
fprintf(fid,'\r\nPSM     : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Psm / t,    objectives.Psm / max_obj,   avgCT.Psm / t);
fprintf(fid,'\r\nGNCCP   : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.Gnccp / t,  objectives.Gnccp / max_obj, avgCT.Gnccp / t);
fprintf(fid, '\r\nBPF_G  : acc: %.2f,  objectives: %.2f, computational time: %.2f', avgAcc.BpfG / t,  objectives.BpfG / max_obj, avgCT.BpfG / t);


fclose(fid);

end