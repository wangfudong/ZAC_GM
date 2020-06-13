function plotMatch_pascal(Map,GT_match,Xnew,Ynew,xinl_ind,yinl_ind,im12,nF,x_shift)
Xn = Xnew(xinl_ind,:);
Yn = Ynew(yinl_ind,:);
LX = length(Xn(:,1));
LY = length(Yn(:,1));
Y_mapped = Map*Yn;
xinl_true = xinl_ind(xinl_ind<=nF);
yinl_true = yinl_ind(yinl_ind<=nF);
xinl_out = (xinl_ind>nF);
yinl_out = (yinl_ind>nF);

hh=figure;set(hh,'position',[500,300,700,430]);
imshow(im12);hold on;
% set(gca,'position',[0.02,0.1,0.96,0.95]);
set(gca,'position',[0.01,0.01,0.98,0.99]);
plot(Xnew(xinl_true,1),Xnew(xinl_true,2),'r.','markersize',20);hold on;
plot(Xn(xinl_out,1),Xn(xinl_out,2),'y*','markersize',8,'linewidth',1.3);hold on;
%text(Xn(:,1),Xn(:,2),num2str(xinl_ind),'color','g','fontsize',13);
plot(Ynew(yinl_true,1)+x_shift,Ynew(yinl_true,2),'r.','markersize',20);hold on;
plot(Yn(yinl_out,1)+x_shift,Yn(yinl_out,2),'y+','markersize',8,'linewidth',2);hold on;
%text(Yn(:,1)+x_shift,Yn(:,2),num2str(yinl_ind),'color','g','fontsize',13);
%plot(Y_mapped(:,1)+x_shift,Y_mapped(:,2),'m+','markersize',15);
%text(Y_mapped(:,1)+x_shift,Y_mapped(:,2),num2str([1:LX]'))

match_pair = matrix2vec(Map);
match_pair1 = [xinl_ind(match_pair(1,:)),yinl_ind(match_pair(2,:))];
for j = 1:length(match_pair(1,:))%LX
%     if j <= nF
        if sum(ismember(GT_match,match_pair1(j,:),'rows'))==1
            line([Xn(match_pair(1,j),1),Yn(match_pair(2,j),1)+x_shift],[Xn(match_pair(1,j),2),Yn(match_pair(2,j),2)],'color','g','linewidth',2);hold on;
        else
            line([Xn(match_pair(1,j),1),Yn(match_pair(2,j),1)+x_shift],[Xn(match_pair(1,j),2),Yn(match_pair(2,j),2)],'color','r','linewidth',2);hold on;
        end
%     else
%         line([Xn(match_pair(1,j),1),Yn(match_pair(2,j),1)+x_shift],[Xn(match_pair(1,j),2),Yn(match_pair(2,j),2)],'linestyle','--','color','m','linewidth',2);hold on;
%     end
        
end

plot(Yn(yinl_out,1)+x_shift,Yn(yinl_out,2),'y+','markersize',8,'linewidth',2);hold on;
plot(Xn(xinl_out,1),Xn(xinl_out,2),'y*','markersize',8,'linewidth',1.3);hold on;
plot(Xnew(xinl_true,1),Xnew(xinl_true,2),'r.','markersize',20);hold on;
plot(Ynew(yinl_true,1)+x_shift,Ynew(yinl_true,2),'r.','markersize',20);hold on;


X_rest = Xnew(setdiff(1:length(Xnew(:,1)),xinl_ind),:);
Y_rest = Ynew(setdiff(1:length(Ynew(:,1)),yinl_ind),:);
plot(X_rest(:,1),X_rest(:,2),'y*','markersize',10,'linewidth',1.5);hold on;
plot(Y_rest(:,1)+x_shift,Y_rest(:,2),'y+','markersize',10,'linewidth',2.5);hold on;