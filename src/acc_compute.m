function [re,pre,xinlrate,yinlrate] = acc_compute(Map,GT,X_inl,Y_inl,nF)

[m,n] = size(Map);
right = 0;
for i = 1:m
    indy = (Map(i,:)==1);
    xmatchy = [X_inl(i),Y_inl(indy)];
    if isempty(Y_inl(indy))
        xmatchy = [X_inl(i),0];
    end
    right = right + (sum(ismember(GT,xmatchy,'rows'))==1);
end
re = right/nF;
pre = right/sum(sum(Map));

xinlrate  = sum(ismember(GT(:,1),X_inl))/length(X_inl);
yinlrate  = sum(ismember(GT(:,2),Y_inl))/length(Y_inl);
