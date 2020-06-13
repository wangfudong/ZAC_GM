function plot_point_reg(X,Y,softmap,inl_idx,inl_idy,X_ind,Y_ind,GT)
X_inl = X(X_ind,:);
Y_inl = Y(Y_ind,:);

plot(X(:,1),X(:,2),'r*',Y(:,1),Y(:,2),'b+','markersize',5);hold on;
plot(X(inl_idx,1),X(inl_idx,2),'r.',Y(inl_idx,1),Y(inl_idy,2),'b.','markersize',20);hold on;


if length(X_ind) <=length(Y_ind)
    %X_mapped1 = softmap*Y_inl;
    %line([X_inl(:,1),X_mapped1(:,1)]',[X_inl(:,2),X_mapped1(:,2)]','color','g');
    cnt_right= 0;
    map_bin = asgHun(softmap);
    X_mapped = map_bin*Y_inl;
    for i = 1:length(X_mapped(:,1))
        [~,c1] = find(map_bin(i,:)==1);
        r1 = X_ind(i);
        c1 = Y_ind(c1);
        if sum(ismember(GT,[r1,c1],'rows'))>0
            line([X_inl(i,1),X_mapped(i,1)]',[X_inl(i,2),X_mapped(i,2)]','color','g');
            cnt_right = cnt_right +1;
        else
            line([X_inl(i,1),X_mapped(i,1)]',[X_inl(i,2),X_mapped(i,2)]','color','r');
        end
    end
    
else
    %Y_mapped1 = softmap'*X_inl;
    %line([Y_inl(:,1),Y_mapped1(:,1)]',[Y_inl(:,2),Y_mapped1(:,2)]','color','g');
    cnt_right= 0;
    map_bin = asgHun(softmap);
    Y_mapped = map_bin'*X_inl;
    for i = 1:length(Y_mapped(:,1))
        [c1,~] = find(map_bin(:,i)==1);
        r1 = Y_ind(i);
        c1 = X_ind(c1);
        if sum(ismember(GT,[c1,r1],'rows'))>0
            line([Y_inl(i,1),Y_mapped(i,1)]',[Y_inl(i,2),Y_mapped(i,2)]','color','g');
            cnt_right = cnt_right +1;
        else
            line([Y_inl(i,1),Y_mapped(i,1)]',[Y_inl(i,2),Y_mapped(i,2)]','color','r');
        end
    end
    
    
end





