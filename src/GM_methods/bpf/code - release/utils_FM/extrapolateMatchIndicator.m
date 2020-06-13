function [ indicator_ext ] = extrapolateMatchIndicator( view, matchList, indicator, search_radius )
% Make an extrapolated GT for evaluation
% output: 1 x nInitialMatches indicator vector
if length(view) == 1
    bSelf = 1;  view1 = 1;  view2 = 1;
else
    bSelf = 0;  view1 = 1;  view2 = 2;
end

if isempty(indicator)
    indicator_ext = ones(1,size(matchList,1));
    fprintf('the given match indicator vector is empty -> all selected...\n');
else
    % find all matches very close to the given matches
    XY_view1= view(view1).feat(matchList(:,1),1:2);
    XY_view2= view(view2).feat(matchList(:,2),1:2); 
    curMatchIdx = find(indicator);
    nCurMatch = length(curMatchIdx);
    indicator_ext = zeros(1,size(matchList,1));
    
    for i=1:nCurMatch
        query_view1 = XY_view1(curMatchIdx(i),1:2);
        query_view2 = XY_view2(curMatchIdx(i),1:2);
        if ~bSelf
            %finds the points within the search radius
            ridx_view1=nearestneighbour(query_view1',XY_view1','Radius',search_radius);
            ridx_view2=nearestneighbour(query_view2',XY_view2','Radius',search_radius);
            closeMatchIdx = intersect(ridx_view1, ridx_view2);
        else
            %finds the points within the search radius
            ridx_view1a=nearestneighbour(query_view1',XY_view1','Radius',search_radius);
            ridx_view2a=nearestneighbour(query_view2',XY_view2','Radius',search_radius);
            ridx_view1b=nearestneighbour(query_view1',XY_view2','Radius',search_radius);
            ridx_view2b=nearestneighbour(query_view2',XY_view1','Radius',search_radius);
            closeMatchIdx = union( intersect(ridx_view1a, ridx_view2a),intersect(ridx_view1b, ridx_view2b));
        end
        indicator_ext(closeMatchIdx) = 1;        
    end
end

% fprintf('input %d matches -> %d extrapolated matches\n', nCurMatch, nnz(indicator_ext) );


