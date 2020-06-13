function plot_interpYX(Map,Xnew,Ynew,xinl_ind,yinl_ind,im12,nF,x_shift)
Xn = Xnew(xinl_ind,:);
Yn = Ynew(yinl_ind,:);
LX = length(Xn(:,1));
LY = length(Yn(:,1));
X_mapped =Map'*Xn;
xinl_true = xinl_ind(xinl_ind<=nF);
yinl_true = yinl_ind(yinl_ind<=nF);
xinl_out = (yinl_ind>nF);
yinl_out = (yinl_ind>nF);
hh=figure;set(hh,'position',[500,300,700,430]);
imshow(im12);hold on;
set(gca,'position',[0.02,0.1,0.96,0.95]);
plot(Xnew(xinl_true,1),Xnew(xinl_true,2),'r.','markersize',15);hold on;
plot(Xn(xinl_out,1),Xn(xinl_out,2),'y*','markersize',15,'linewidth',2);hold on;
text(Xn(:,1),Xn(:,2),num2str(xinl_ind),'color','g','fontsize',13);
plot(Ynew(yinl_true,1)+x_shift,Ynew(yinl_true,2),'r.','markersize',15);hold on;
plot(Yn(yinl_out,1)+x_shift,Yn(yinl_out,2),'y+','markersize',15,'linewidth',2);hold on;
text(Yn(:,1)+x_shift,Yn(:,2),num2str(yinl_ind),'color','g','fontsize',13);
plot(X_mapped(:,1),X_mapped(:,2),'m+','markersize',15);
%text(X_mapped(:,1),X_mapped(:,2),num2str([1:LY]'));
for j = 1:LY
    line([Yn(j,1)+x_shift,X_mapped(j,1)],[Yn(j,2),X_mapped(j,2)],'color',rand(1,3),'linewidth',2);
end