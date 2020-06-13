function [X_inl,Y_inl] = out_re(X_index,Y_index,Map,th,mass,type)
% outlier identification and removal

if nargin == 5 || type ==1
    [LX,LY] = size(Map);
    [row_sum,row_ind] = sort(sum(Map,2),'descend');
    [col_sum,col_ind] = sort(sum(Map),'descend');
    col_sum = col_sum';col_ind = col_ind';
    
    if LX < LY
        row_sum(LX+1:LY,:) = 0;
        row_ind(LX+1:LY,:) = (LX+1):1:LY;
    elseif LX > LY
        col_sum(LY+1:LX,:) = 0;
        col_ind(LY+1:LX,:) = (LY+1):1:LX;
    end
    
    [INDEX, C] = kmeans([row_sum,col_sum], 2);
    CC = sum(C,2);

% %     tic;
%     [C, INDEX] = vl_kmeans([row_sum,col_sum]', 2);% if vlfeat-0.9.20-bin is used 
%     CC = sum(C,1);
% %     toc;

    if CC(1) > CC(2)
        INDEX_out = find(INDEX==2);
        %INDEX_inl = find(INDEX==1);
    elseif CC(1) < CC(2)
        INDEX_out = find(INDEX==1);
        %INDEX_inl = find(INDEX==2);
    elseif CC(1) == CC(2)
        M = sum(INDEX==1);
        L = sum(INDEX==2);
        if M >= L
            INDEX_out = find(INDEX==1);
        else
            INDEX_out = find(INDEX==2);
        end
        
    end
    
    if length(INDEX_out) > (max(LX,LY)-mass)
        INDEX_out = INDEX_out((length(INDEX_out)-max(LX,LY)+mass+1):end);
        %     row_out = (row_sum(INDEX_out) <= th);
        %     col_out = (col_sum(INDEX_out) <= th);
        %     row_out = INDEX_out(row_out);
        %     col_out = INDEX_out(col_out);
        % else
        %     row_out = INDEX_out;
        %     col_out = INDEX_out;
    end
    
    row_out = (row_sum(INDEX_out) <= th);
    col_out = (col_sum(INDEX_out) <= th);
    row_out = INDEX_out(row_out);
    col_out = INDEX_out(col_out);
    
    row_out = row_ind(row_out,:);
    col_out = col_ind(col_out,:);
    
%     row_inl = setdiff(row_ind,row_out);
%     col_inl = setdiff(col_ind,col_out);
    row_inl = sort(row_ind(~ismember(row_ind,row_out)));
    col_inl = sort(col_ind(~ismember(col_ind,col_out)));

    X_inl = X_index(row_inl);
    Y_inl = Y_index(col_inl);
    
%     X_out = setdiff(X_index,X_inl);
%     Y_out = setdiff(Y_index,Y_inl);
    
else
    
    r1 = sum(Map,2);
    c1 = sum(Map,1);
    rc_map = r1*c1;
    map = kLAP_Hun(-rc_map,mass);
    [r_ind,c_ind] = find(map==1);    
    
    
    X_inl = X_index(r_ind);
    Y_inl = Y_index(c_ind);
    X_inl = sort(X_inl);
    Y_inl = sort(Y_inl);
    
%     X_out = setdiff(X_index,X_inl);
%     Y_out = setdiff(Y_index,Y_inl);
    
end













