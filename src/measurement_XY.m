function reg_err = measurement_XY(X_new,Y_tmp,X_tmp,cor_GT)
[~,max_boundx] = normalize_point(X_tmp(cor_GT(:,1),:),1);
[~,max_boundy] = normalize_point(Y_tmp(cor_GT(:,2),:),1);
max_bound = max(max_boundx,max_boundy);

reg_err = mean(sqrt(sum((X_new(cor_GT(:,1),:) - Y_tmp(cor_GT(:,2),:)).^2,2)),1)/max_bound;
