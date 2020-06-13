function Map_all = Map_by_TopK(map_soft,map_bin,topk)
[m,n] = size(map_soft);
map_scores = map_soft.*map_bin;
[row,col] = find(map_bin==1);
max_scores = zeros(length(row),1);
for ii = 1:length(row)
    max_scores(ii) = map_scores(row(ii),col(ii));
end
rand_order = randperm(length(row));
max_scores = max_scores(rand_order);
[~,inds] = sort(max_scores,'descend');
inds_max = rand_order(inds(1:topk));
Map_all = zeros(m,n);
for j = 1:length(inds_max)
    Map_all(row(inds_max(j)),col(inds_max(j)))=1;
end