function para = RRWM_nonrigid_2d(X,Y,opt_reg,opt_gm)

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

xinl = cor_GT(:,1);
yinl = cor_GT(:,2);
xout = setdiff(1:m,xinl)';
yout = setdiff(1:n,yinl)';

Ws = zeros(2,m,2);

reg_err = measurement_XY(X,Y,X,opt_gm.GT);
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

Ct = ones(m,n);
[ind, ind_m] = find(Ct);
group1 = zeros(size(ind, 1), m);
group2 = zeros(size(ind, 1), n);
for i = 1:size(ind, 1)
    group1(i, ind(i)) = 1;
    group2(i, ind_m(i)) = 1;
end
group1 = logical(group1);
group2 = logical(group2);

GXX = nonrigid_ker(X_tmp,opt_reg.nonrigid_sig,'rbf');

while reg_cnt <= max_reg && reg_diff >= 1.0e-5 %%regist_err >= 1.0e-5 
    
    opt_reg.SC_or_DIS = unary_term(reg_cnt);
    
    Pts{1,1} = X';
    Pts{1,2} = Y';
    parKnl = st('alg', 'pas'); % type of affinity: only edge distance
    parG = st('link', opt_gm.full_or_del); % Delaunay triangulation for computing the graphs
    gphs = newGphUs(Pts, parG);
    
    [~, KQ] = conKnlGphPQU(gphs, parKnl);
    
    if opt_reg.SC_or_DIS == 0
        KP =  M_shape(Pts{1,1}',Pts{1,2}',1/8,2,opt_reg.rota);
        KP = exp(-KP/0.1);
    else
        KP = M_points(Pts{1,1}',Pts{1,2}');
        KP = exp(-KP/0.1);
    end
    
    K = conKnlGphKU(KP, KQ, gphs);
    
    RRWM_soft = RRWM(K, group1, group2);
    RRWM_map = greedyMapping(RRWM_soft, group1, group2);
    RRWM_map = reshape(RRWM_map, m, n);
    
    RRWM_scores = reshape(RRWM_soft, m, n).*RRWM_map;
    [row,col] = find(RRWM_map==1);
    max_scores = zeros(length(row),1);
    for ii = 1:length(row)
        max_scores(ii) = RRWM_scores(row(ii),col(ii));
    end
    rand_order = randperm(length(row));
    max_scores = max_scores(rand_order);
    [~,inds] = sort(max_scores,'descend');
    inds_max = rand_order(inds(1:opt_gm.topk));
    Map_all = zeros(size(RRWM_map));
    for j = 1:length(inds_max)
        Map_all(row(inds_max(j)),col(inds_max(j)))=1;
    end
    
    NP = sum(Map_all(:));
    
    sig2 = max(sum(sum(Map_all.*(bsxfun(@minus,X(:,1),Y(:,1)').^2 + bsxfun(@minus,X(:,2),Y(:,2)').^2)))/(2*NP),1.0e-10);
    W1 = ((opt_reg.lambda*sig2*eye(size(GXX))+diag(sum(Map_all,2))*GXX))\(Map_all*Y-diag(sum(Map_all,2))*X_tmp);
    
    Ws(reg_cnt,:,:) = W1;
    
    X_old = X;
    X = nonrigid_kernel_trans(X_tmp,W1,GXX,eye(2,2),[]);
    
    reg_err = measurement_XY(X,Y_tmp,Y_tmp,opt_gm.GT);
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
para.W = Ws;
para.reg_Err = reg_errs;
para.reg_Diff = reg_diffs;
para.X_new = X;
para.X_tmp = X_tmp;
para.Y_tmp = Y_tmp;

para.ux = ux;
para.uy = uy;
