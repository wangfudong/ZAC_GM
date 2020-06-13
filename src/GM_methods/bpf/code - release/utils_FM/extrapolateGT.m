function [ indicator_GT_ext ] = extrapolateGT( view, matchList, GTList, thresh_dist_GT )
% Make an extrapolated GT for evaluation
% output: 1 x nInitialMatches indicator vector

if length(view) == 1
    bSelf = 1;  view1 = 1;  view2 = 1;
else
    bSelf = 0;  view1 = 1;  view2 = 2;
end

if isempty(GTList)
    indicator_GT_ext = ones(1,size(matchList,1));
    fprintf('Ground truth not included -> all matches considered as True... \n');
else
    % find all matches close to GTs among the initial matches
    XY_view1= view(view1).feat(matchList(:,1),1:2);
    XY_view2= view(view2).feat(matchList(:,2),1:2); 
    %matchList_GT = matchList(find(GT),:);
    nGTList = size(GTList,1);
    indicator_GT_ext = zeros(1,size(matchList,1));
    
    for i=1:nGTList
        query_view1 = view(view1).feat(GTList(i,1),1:2);
        query_view2 = view(view2).feat(GTList(i,2),1:2); 
        if ~bSelf
            %finds the points within the search radius
            ridx_view1=nearestneighbour(query_view1',XY_view1','Radius',thresh_dist_GT);
            ridx_view2=nearestneighbour(query_view2',XY_view2','Radius',thresh_dist_GT);
            trueIdx = intersect(ridx_view1, ridx_view2);
        else
            %finds the points within the search radius
            ridx_view1a=BruteSearchMex(query_view1',XY_view1','Radius',thresh_dist_GT);
            ridx_view2a=BruteSearchMex(query_view2',XY_view2','Radius',thresh_dist_GT);
            ridx_view1b=BruteSearchMex(query_view1',XY_view2','Radius',thresh_dist_GT);
            ridx_view2b=BruteSearchMex(query_view2',XY_view1','Radius',thresh_dist_GT);
            trueIdx = union( intersect(ridx_view1a, ridx_view2a),intersect(ridx_view1b, ridx_view2b));
        end
        indicator_GT_ext(trueIdx) = 1;        
    end
end



