function cdata = conRrwmGphs(dataPath, cImg, outliers)

%% load data
datafile = sprintf('%s/matchData/fi_%02da+%02db.mat', dataPath, cImg, cImg);
load(datafile);
    
%% affinity matrix
affinity_max = 50;               % maximum value of affinity 
cdata.affinityMatrix = max(affinity_max - cdata.distanceMatrix,0); % dissimilarity -> similarity conversion
cdata.affinityMatrix(1:(length(cdata.affinityMatrix)+1):end) = 0; % diagonal zeros

curMatchList = cell2mat({cdata.matchInfo(:).match }');
num1 = size(cdata.view(1).feat, 1);     % number of feature points
num2 = size(cdata.view(2).feat, 1);     % number of feature points


ind = 1;
for nOut = outliers

    %% select randomly feature points
    ptList1 = 1:num1;
    ptList2 = 1:num2;

    % ground-truth points
    ptSel1 = unique(cdata.GT(:,1));
    ptSel2 = unique(cdata.GT(:,2));
    nGt1 = size(ptSel1, 1);
    nGt2 = size(ptSel2, 1);

    % outliers
    if nOut > 0
        ptList1(ptSel1) = 0;
        ptList2(ptSel2) = 0;
        ptList1 = ptList1(ptList1 > 0);
        ptList2 = ptList2(ptList2 > 0);

        nOut1 = min(nOut, size(ptList1, 2));
        nOut2 = min(nOut, size(ptList2, 2));
        ptSel1 = [ptSel1; ptList1(randperm(num1 - nGt1, nOut1))'];
        ptSel2 = [ptSel2; ptList2(randperm(num2 - nGt2, nOut2))'];
    end
    
    %% filter unneccesary matches
    bPointSel1 = zeros(num1, 1);
    bPointSel1(ptSel1) = 1;
    bPointSel1 = logical(bPointSel1);
    
    bPointSel2 = zeros(num2, 1);
    bPointSel2(ptSel2) = 1;
    bPointSel2 = logical(bPointSel2);

    bMatchSel = zeros(size(curMatchList, 1), 1);
    for i = 1:size(curMatchList, 1)
        if bPointSel1(curMatchList(i,1)) && bPointSel2(curMatchList(i,2))
            bMatchSel(i) = 1;
        end
    end
    bMatchSel = logical(bMatchSel);
    
    selMatchList = curMatchList(bMatchSel, :);
    K    = cdata.affinityMatrix(bMatchSel, bMatchSel);
    X_GT = zeros(size(selMatchList, 1), 1);
    for i = 1:size(selMatchList, 1)
        for k = 1:size(cdata.GT, 1)
            if selMatchList(i,1) == cdata.GT(k,1) && selMatchList(i,2) == cdata.GT(k,2)
                X_GT(i) = 1;
                break;
            end
        end
    end
    
    group1 = cdata.group1(bMatchSel, :);
    group1 = group1(:,any(group1,1));

    group2 = cdata.group2(bMatchSel, :);
    group2 = group2(:,any(group2,1));
    
    cdata.GRAPH{ind}.ptSel1   = ptSel1;
    cdata.GRAPH{ind}.ptSel2   = ptSel2;
    cdata.GRAPH{ind}.selMatchList = selMatchList;
    cdata.GRAPH{ind}.X_GT     = X_GT;
    cdata.GRAPH{ind}.K        = K;
    cdata.GRAPH{ind}.group1   = group1;
    cdata.GRAPH{ind}.group2   = group2;
    
    ind = ind + 1;
end

end