function [Map,output] = ZAC_vgg_r(X,Y,opt)

DX = M_points(X,X);
DY = M_points(Y,Y);

max_fact = max([max(DX(:)),max(DY(:))]);
DX = DX/max_fact;
DY = DY/max_fact;


SX = 1./(DX+eye(size(DX)))-eye(size(DX));SX = SX/max(SX(:));
SY = 1./(DY+eye(size(DY)))-eye(size(DY));SY = SY/max(SY(:));


if strcmp(opt.del_or_full,'del') == 1
    sx = nei_del(X);
    sx = sx + sx';
    sy = nei_del(Y);
    sy = sy + sy';
    
    DX_tmp = DX.*sx;
    DY_tmp = DY.*sy;
    SX_tmp = SX.*sx;
    SY_tmp = SY.*sy;
    
else
    DX_tmp = DX;
    DY_tmp = DY;
    SX_tmp = SX;
    SY_tmp = SY;
end

LX = length(X(:,1));
LY = length(Y(:,1));
X_inl_ind = (1:LX)';
Y_inl_ind = (1:LY)';
% XX_inl = X;
% YY_inl = Y;
descx = opt.descx;
descy = opt.descy;

cnt = 1;
stop = 1;

while cnt <= opt.remove_time && stop > 0
   
    DE1 = descx(:,X_inl_ind);
    DE2 = descy(:,Y_inl_ind);
    %M = hist_cost_2(DE1',DE2');M = M/max(M(:));
    M = Dis_sift(DE1',DE2');

    DX = DX_tmp(X_inl_ind,:); DX = DX(:,X_inl_ind);
    DY = DY_tmp(Y_inl_ind,:); DY = DY(:,Y_inl_ind);

    SX = SX_tmp(X_inl_ind,:); SX = SX(:,X_inl_ind);
    SY = SY_tmp(Y_inl_ind,:); SY = SY(:,Y_inl_ind);
       
    
    sigx = max(std(DX(DX(:)>0))^2,0.01);
    sigy = max(std(DY(DY(:)>0))^2,0.01);
    sig = [sigx,sigy];
   % sig = [0.1,0.1].^2;
    A = exp(-DX.^2/sig(1)).*(double(DX>0) + eye(size(DX)));
    B = exp(-DY.^2/sig(2)).*(double(DY>0) + eye(size(DY)));

    
    Map_ini = ones(LX,LY)/max(LX,LY);
    
    [Map,fw_output] = FW_minimu_r(Map_ini,M,SX,SY,A,B,opt);
    if opt.r_fixed > 0
        opt.mass_r = sum(Map(:));
    end
    mass1 = min(max(round(sum(Map(:))),5),min(LX,LY));
    
    if opt.output > 0
        output.optim_output(cnt).optim_output = fw_output;
        output.Maps(cnt).Maps = Map;
        output.mass_r(cnt) = sum(sum(Map));
    end
    
    if opt.remove_time > 1 && cnt < opt.remove_time
        if cnt <= 1
            th = 0.2;
        else
            th = 0.2;
        end
        
        X_inl_tmp = X_inl_ind;
        Y_inl_tmp = Y_inl_ind;
        
        if ~isfield(opt,'out_re_type')
            [X_inl_ind,Y_inl_ind] = out_re(X_inl_ind,Y_inl_ind,Map,th,mass1);
        else
            [X_inl_ind,Y_inl_ind] = out_re(X_inl_ind,Y_inl_ind,Map,th,mass1,opt.out_re_type);
        end
        
%         XX_inl = X(X_inl_ind,:);
%         YY_inl = Y(Y_inl_ind,:);
        LX = length(X_inl_ind);
        LY = length(Y_inl_ind);
        
        X_stop = sum(~ismember(X_inl_tmp,X_inl_ind));
        Y_stop = sum(~ismember(Y_inl_tmp,Y_inl_ind));
        stop = (X_stop+Y_stop);
        cnt = cnt+1;
        
        if opt.output > 0
            output.inl(cnt).xinl = X_inl_ind;
            output.inl(cnt).yinl = Y_inl_ind;
        end
        
    else 
        stop = 0;
    end
    
    output.X_inl = X_inl_ind;
    output.Y_inl = Y_inl_ind;
    
end




