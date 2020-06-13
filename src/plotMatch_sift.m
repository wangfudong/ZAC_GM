function [dis] = plotMatch_sift(Map_bin,Xnew,Ynew,xinl_ind,yinl_ind,im12,x_shift,H,pix_tol)
Xn = Xnew(xinl_ind,:);
Yn = Ynew(yinl_ind,:);

match_pairs = matrix2vec(Map_bin);
dis = 0;

hh=figure;set(hh,'position',[500,300,700,430]);
imshow(im12);hold on;
set(gca,'position',[0.005,0.005,0.99,0.995]);
plot(Xnew(:,1),Xnew(:,2),'y+','markersize',10,'linewidth',3);hold on;
plot(Xn(:,1),Xn(:,2),'r.','markersize',20);hold on;
plot(Ynew(:,1)+x_shift,Ynew(:,2),'y+','markersize',10,'linewidth',3);hold on;
plot(Yn(:,1)+x_shift,Yn(:,2),'r.','markersize',20);hold on;

XM = Xn(match_pairs(1,:),:);
X_matched = Yn(match_pairs(2,:),:);
XM_mapped = zeros(3,length(XM(:,1)));
ttt = zeros(1,length(XM(:,1)));

for i = 1:length(XM(:,1))
    XM_mapped(:,i) = H*[XM(i,:),1]';
    XM_mapped(:,i) = XM_mapped(:,i)/XM_mapped(3,i);
    Dis = sqrt((XM_mapped(1,i) - X_matched(i,1)).^2 + (XM_mapped(2,i) - X_matched(i,2)).^2);
    dis(i) = Dis;
    if Dis <= pix_tol
        line([XM(i,1),X_matched(i,1)+x_shift],[XM(i,2),X_matched(i,2)],'color','g','linewidth',2);hold on;
    else

        line([XM(i,1),X_matched(i,1)+x_shift],[XM(i,2),X_matched(i,2)],'color','r','linewidth',2);hold on;
    end
    ttt(i) = (Dis <= pix_tol);
end
