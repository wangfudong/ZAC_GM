function [X_inl,X_out] = pascal_pick_outlier(X_all,inl_index,out_index,out_num)
X_inl = X_all(inl_index,:);
if isempty(out_index) == 0
    X_out = X_all(out_index,:);
else
    X_rest = X_all(length(X_inl(:,1))+1:end,:);
    out_inds = randperm(length(X_all(:,1))-length(X_inl(:,1)));% add outliers randomly
    out_inds = out_inds(1:min(out_num,length(out_inds)));
    X_out = X_rest(out_inds,:);
end