function para = FRGM_similar_2d_fast(X,Y,opt_reg,opt_frgm)

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

if isfield(opt_frgm,'GT')
    cor_GT = opt_frgm.GT;
else
    cor_GT = [1:min(m,n);1:min(m,n)]';
end

max_reg = opt_reg.reg_maxiter;
reg_errs = zeros(1,1);
reg_diffs = zeros(1,1);

xinl = cor_GT(:,1);
yinl = cor_GT(:,2);
xout = setdiff(1:m,xinl)';
yout = setdiff(1:n,yinl)';

Rs = zeros(max_reg,dim1,dim1);
Ts = zeros(max_reg,2);
Ss = zeros(max_reg,1);

reg_err = measurement_XY(X,Y,X,opt_frgm.GT);
reg_diff = reg_err;
reg_cnt = 1;
unary_term = opt_reg.unary_terms;

file_write = 'save/';
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

GT = opt_frgm.GT;
asgT.X = zeros(m,n);
for i = 1:length(GT(:,1))
    asgT.X(GT(i,1),GT(i,2)) = 1;
end

while reg_cnt <= max_reg && reg_diff >= 1.0e-4 %%regist_err >= 1.0e-5 
    
    opt_reg.SC_or_DIS = unary_term(reg_cnt);
    if  opt_reg.SC_or_DIS == 0
        [M,~,~] = M_shape(X,Y,1/8,2,opt_frgm.rota);
        opt_frgm.geofunc = '1.21';
    else
        M = M_points(X,Y);
        M = M/max(M(:));
        opt_frgm.geofunc = '1.21';'0.11';
    end
    opt_frgm.convex = unary_term(reg_cnt);
    
    Map_ini = asgHun(-M);
    DXX = M_points(X,X);
    SX = 1./(DXX+eye(size(DXX)))-eye(size(DXX));
    SX = SX/max(SX(:));
    FRGM_scores = FW_Min(X,Y,M,SX,Map_ini,opt_frgm,1);
    
    [row,col] = find(asgHun(FRGM_scores)==1);
    max_scores = zeros(length(row),1);
    for ii = 1:length(row)
        max_scores(ii) = FRGM_scores(row(ii),col(ii));
    end
    rand_order = randperm(length(row));
    max_scores = max_scores(rand_order);
    [~,inds] = sort(max_scores,'descend');
    inds_max = rand_order(inds(1:opt_frgm.topk));
    Map_all = zeros(size(FRGM_scores));
    for j = 1:length(inds_max)
        Map_all(row(inds_max(j)),col(inds_max(j)))=1;
    end
  
    [R1,s1,t1] = rigid_parameter(X,Y,Map_all);
        
    Rs(reg_cnt,:,:) = R1;
    Ts(reg_cnt,:) = t1;
    Ss(reg_cnt) = s1;
    
    X_old = X;
    X = s1*X*R1 + repmat(t1,m,1);
    
    reg_err = measurement_XY(X,Y_tmp,Y_tmp,opt_frgm.GT);
    reg_diff = measurement_XY(X_old,X,X_old,[1:m;1:m]');
    reg_errs(reg_cnt) = reg_err;
    reg_diffs(reg_cnt) = reg_diff;

    
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
para.R = Rs;
para.T = Ts;
para.S = Ss;
para.reg_Err = reg_errs;
para.reg_Diff = reg_diffs;
para.X_new = X;
para.X_tmp = X_tmp;
para.Y_tmp = Y_tmp;


