function [X2Y_correct,Y2X_correct,xre,xpre,yre,ypre] = match_result_sift(X,Y,X_inlnum,Y_inlnum,match_pairs,H,tol)
% X:m*2
% Y:n*2
% match_pairs: 2*k
% H:3*3

X2Y_correct = 0;
XM = X(match_pairs(1,:),:);
X_matched = Y(match_pairs(2,:),:);
XM_mapped = zeros(3,length(XM(:,1)));
for i = 1:length(XM(:,1))
    XM_mapped(:,i) = H*[XM(i,:),1]';
    XM_mapped(:,i) = XM_mapped(:,i)/XM_mapped(3,i);
    Dis = sqrt((XM_mapped(1,i) - X_matched(i,1)).^2 + (XM_mapped(2,i) - X_matched(i,2)).^2);
    if Dis <= tol
        X2Y_correct = X2Y_correct + 1;
    end
end


Y2X_correct = 0;
YM = Y(match_pairs(2,:),:);
Y_matched = X(match_pairs(1,:),:);
YM_mapped = zeros(3,length(YM(:,1)));
for i = 1:length(YM(:,1))
    YM_mapped(:,i) = H\[YM(i,:),1]';
    YM_mapped(:,i) = YM_mapped(:,i)/YM_mapped(3,i);
    Dis = sqrt((YM_mapped(1,i) - Y_matched(i,1)).^2 + (YM_mapped(2,i) - Y_matched(i,2)).^2);
    if Dis <= tol
        Y2X_correct = Y2X_correct + 1;
    end
end

xre = X2Y_correct/X_inlnum;
xpre = X2Y_correct/length(match_pairs(1,:));

yre = Y2X_correct/Y_inlnum;
ypre = Y2X_correct/length(match_pairs(1,:));




