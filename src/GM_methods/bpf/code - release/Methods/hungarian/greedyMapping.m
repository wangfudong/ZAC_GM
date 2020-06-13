function Xd = greedyMapping(score, group1, group2)
% greedy discretization using the group1 and group2
% group1 = # of matches by ( # of groups in feat1 )
% group2 = # of matches by ( # of groups in feat2 )

Xd = zeros(length(score),1); % discretized solution (binary vector)
[ max_value max_ind ] = max(score);
while max_value > 0
    Xd(max_ind) = 1;
    group1_idx = find(group1(max_ind,:));
    score(find(group1(:,group1_idx))) = 0;
    group2_idx = find(group2(max_ind,:));
    score(find(group2(:,group2_idx))) = 0;
    [ max_value max_ind ] = max(score);
end
