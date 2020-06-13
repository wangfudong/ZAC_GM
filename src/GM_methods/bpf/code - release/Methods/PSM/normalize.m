function p = normalize(q, group)
% normalized by group1
[~, group_id] = find(group);
nGroups = max(group_id);
for gid = 1:nGroups
    [match_id, ~] = find(group(:,gid));
    norm_q = sum(q(match_id)) + eps;
    q(match_id) = q(match_id) / norm_q;
end
p = q;
end