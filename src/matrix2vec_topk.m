function M_topks = matrix2vec_topk(softmap,rates)

[m,n] = size(softmap);
M_topks = zeros(length(rates),m,n);

if m <= n
    max_scores = zeros(m,1);
    max_scores_id = zeros(m,2);
    for i = 1:m
        max_scores(i) = max(softmap(i,:));
        if max_scores(i) > 0
            ci = find(softmap(i,:)==max_scores(i));
            max_scores_id(i,:) = [i,ci(1)];
        end
    end
else
    max_scores = zeros(n,1);
    max_scores_id = zeros(n,2);
    for i = 1:n
        max_scores(i) = max(softmap(:,i));
        if max_scores(i) > 0
            ci = find(softmap(:,i)==max_scores(i));
            max_scores_id(i,:) = [ci(1),i];
        end
    end
        
end

[~,inds] = sort(max_scores,'descend');
max_scores_id1 = max_scores_id(inds,:);

m1 = length(rates);

for i = 1:m1
    TopK = round(rates(i)*min(m,n));
    if TopK==0
        error('k is too small!');
    end
    for j = 1:TopK
        if max_scores_id1(j,1) > 0 && max_scores_id1(j,2) >0
            M_topks(i,max_scores_id1(j,1),max_scores_id1(j,2)) = 1;
        end
    end
    
end

% 
% Map = zeros(2,length(M(:,1)));
% for i = 1:length(M(:,1))
%     indi = find(M(i,:)==1);
%     if isempty(indi)~=1
%         Map(1,i) = i;
%         Map(2,i) = indi;
%     end
% end
% zero_ind = ismember(Map',[0,0],'rows');
% match_pair = Map(:,zero_ind==0);