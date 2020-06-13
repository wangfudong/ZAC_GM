function [PRs,max_scores,max_scores_id] = TopK_PR(softmap,gt,nF,rates)

[m,n] = size(softmap);

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
PRs = zeros(2,m1);
[r,c]=find(gt==1);

for i = 1:m1
    TopK = round(rates(i)*min(m,n));
    if TopK==0
        error('k is too small!');
    end
    right = 0;
    for j=1:TopK
    right = right+ sum(ismember([r,c],max_scores_id1(j,:),'rows'));
    end
    PRs(1,i) = right/nF;
    PRs(2,i) = right/TopK;
    
end


