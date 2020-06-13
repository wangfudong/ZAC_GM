function para = ZAC_nonrigid_2d_r(X,Y,opt_reg,opt_gm)

[m,dim1] = size(X);
[n,dim2] = size(Y);

if dim1~=2 || dim1~=dim2
    error('Only for 2d case!');
end

if opt_reg.reg_normalize > 0
    [X,sx0,ux] = normalize_point(X,1);
    [Y,sy0,uy] = normalize_point(Y,1);
else
    sx0 = 1;ux = [0,0];
    sy0 = 1;uy = [0,0];
end

X_tmp = X;
Y_tmp = Y;

if isfield(opt_gm,'GT')
    cor_GT = opt_gm.GT;
else
    cor_GT = [1:min(m,n);1:min(m,n)]';
end

max_reg = opt_reg.reg_maxiter;
reg_errs = zeros(1,1);
reg_diffs = zeros(1,1);

GM_OUT = struct([]);
for o = 1:max_reg
    GM_OUT(o).output_gm = [];
end

xinl = cor_GT(:,1);
yinl = cor_GT(:,2);
xout = setdiff(1:m,xinl)';
yout = setdiff(1:n,yinl)';

nF_out = zeros(1,2);
Map_out = zeros(2,m,n);
Ws = zeros(2,m,2);
Sig2 = zeros(2,1);

reg_err = measurement_XY(X,Y,X,opt_gm.GT);
reg_diff = reg_err;
reg_cnt = 1;
unary_term = opt_reg.unary_terms;

file_write = 'E:\computer_vision\code\point-reg\Test\similar\';
if opt_reg.display > 0
    hh = figure('color','white');subplot('position',[0.02,0.05,0.45,0.9]);
    set(hh,'position',[500, 450 ,700 ,300]);
    plot(X(xinl,1),X(xinl,2),'r.',Y(yinl,1),Y(yinl,2),'b.','markersize',20);hold on;
    plot(X(xout,1),X(xout,2),'r+',Y(yout,1),Y(yout,2),'b+');
    hold off;
    axis off;
    % axis equal;
    subplot('position',[0.52,0.05,0.45,0.9]);   
    % set(gca,'position',[0,1,0.1,0.8,0.9]);
    plot(X(xinl,1),X(xinl,2),'r.',Y(yinl,1),Y(yinl,2),'b.','markersize',20);hold on;
    plot(X(xout,1),X(xout,2),'r+',Y(yout,1),Y(yout,2),'b+');
    hold off;
    axis off;
    % htext = text(-0.9,0.5,['Registration time = 0'],'fontsize',15);
    drawnow;
    
    if opt_reg.write > 0
        frame = getframe(hh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        %imwrite(imind,cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
        %set(htext,'visible','off');
        imwrite(imind,cm,[file_write '56.jpg'],'jpg');
    end
end

GXX = nonrigid_ker(X_tmp,opt_reg.nonrigid_sig,'rbf');

while reg_cnt <= 5 || (reg_cnt <= max_reg && reg_diff >= 1.0e-4) %%regist_err >= 1.0e-5 
    
    nF = nF_by_kmeans(X,Y,unary_term(reg_cnt),0);
%     nF = max(nF,10);
%     opt_gm.mass = nF;
    nF_out(reg_cnt) = nF;
    
    opt_reg.SC_or_DIS = unary_term(reg_cnt);
    [Map,output_gm] = ZAC_reg_r(X,Y,opt_gm,opt_reg);
    NP = sum(Map(:));
    
    Map_all = Map_recover(Map,m,n,output_gm.X_inl,output_gm.Y_inl);
   
    sig2 = max(sum(sum(Map_all.*(bsxfun(@minus,X(:,1),Y(:,1)').^2 + bsxfun(@minus,X(:,2),Y(:,2)').^2)))/(2*NP),1.0e-10);
    W1 = ((opt_reg.lambda*sig2*eye(size(GXX))+diag(sum(Map_all,2))*GXX))\(Map_all*Y-diag(sum(Map_all,2))*X_tmp);
    
    X_old = X;
    X = nonrigid_kernel_trans(X_tmp,W1,GXX,eye(2,2),[]);
    
    reg_err = measurement_XY(X,Y_tmp,Y_tmp,opt_gm.GT);
    reg_diff = measurement_XY(X_old,X,X_old,[1:m;1:m]');
    reg_errs(reg_cnt) = reg_err;
    reg_diffs(reg_cnt) = reg_diff;
    
    if opt_gm.output > 0
        GM_OUT(reg_cnt).output_gm = output_gm;
        Map_out(reg_cnt,:,:) = Map_all;
        Ws(reg_cnt,:,:) = W1;
        Sig2(reg_cnt) = sig2;
    end
    
   disp(['reg_err at iter ' num2str(reg_cnt) ': ' num2str(reg_err)]);
   disp(['reg_diff at iter ' num2str(reg_cnt) ': ' num2str(reg_diff)]);
    
    reg_cnt = reg_cnt + 1;
    
    if opt_reg.display > 0
        plot(Y(yinl,1),Y(yinl,2),'b.','markersize',10); hold on;
        plot(Y(yout,1),Y(yout,2),'b+'); hold on;
        plot(X(xinl,1),X(xinl,2),'ro','linewidth',1.5,'markersize',6);hold on;
        plot(X(xout,1),X(xout,2),'r+');
        % hold on,plot(Y(LX+1:end,1),Y(LX+1:end,2),'g.','markersize',10);
        %axis equal;
        hold off;
        axis off;
        %htext = text(-0.9,0.55,['Registration = ' num2str(reg_cnt)],'fontsize',15);
        drawnow;
        
        if opt_reg.write > 0
            
            frame = getframe(hh);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            %imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.4);
            % set(htext,'visible','off');
            
            imwrite(imind,cm,[file_write num2str(reg_cnt+55) '.jpg'],'jpg');
        end
        
    end
    
end

reg_cnt = reg_cnt -1;
disp(['registration_err  at iter ' num2str(reg_cnt) ': ' num2str(reg_err)]);
disp(['registration_diff at iter ' num2str(reg_cnt) ': ' num2str(reg_diff)]);

para.sx0 = sx0;
para.sy0 = sy0;
para.ux = ux;
para.uy = uy;
if opt_gm.output > 0
para.GM_OUT = GM_OUT;
para.Maps = Map_out;
para.W = Ws;
end
para.Sig2 = Sig2;
para.reg_Err = reg_errs;
para.reg_Diff = reg_diffs;
para.X_new = X;
para.X_tmp = X_tmp;
para.Y_tmp = Y_tmp;
para.nF_out = nF_out;


