function [Map,output] = ZAC_reg_r(X,Y,opt_gm,opt_reg)

DX = M_points(X,X);
DY = M_points(Y,Y);

max_fact = max([max(DX(:)),max(DY(:))]);
DX = DX/max_fact;
DY = DY/max_fact;

SX = 1./(DX+eye(size(DX)))-eye(size(DX));SX = SX/max(SX(:));
SY = 1./(DY+eye(size(DY)))-eye(size(DY));SY = SY/max(SY(:));

if strcmp(opt_gm.del_or_full,'del') == 1
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
XX_inl = X;
YY_inl = Y;

if opt_gm.output >0
    output.optim_output = struct([]);
    output.Maps = struct([]);
    output.inl = struct([]);
    output.inl(1).xinl = X_inl_ind;
    output.inl(1).yinl = Y_inl_ind;
    output.mass_r = [];
end

cnt = 1;
stop = 1;

while cnt <= opt_gm.remove_time && stop > 0
    
    if opt_reg.SC_or_DIS == 0
    M = M_shape(XX_inl,YY_inl,1/8,2,opt_reg.rota);
    else
        M = M_points(XX_inl,YY_inl);
        M = M/max(M(:));
    end
    DX = DX_tmp(X_inl_ind,:); DX = DX(:,X_inl_ind);
    DY = DY_tmp(Y_inl_ind,:); DY = DY(:,Y_inl_ind);
    
    sigx = std(DX(DX(:)>0))^2;
    sigy = std(DY(DY(:)>0))^2;
    sig = 1*[sigx,sigy];

    A = exp(-DX.^2/sig(1)).*(double(DX>0) + eye(size(DX)));
    B = exp(-DY.^2/sig(2)).*(double(DY>0) + eye(size(DY)));
    if strcmp(opt_reg.type,'nonrigid')==1
%         SX = double_stochastic(A-eye(size(A)),1000);
%         SY = double_stochastic(B-eye(size(B)),1000);
        SX = SX_tmp(X_inl_ind,:); SX = SX(:,X_inl_ind);
        SY = SY_tmp(Y_inl_ind,:); SY = SY(:,Y_inl_ind);
    else
        SX = SX_tmp(X_inl_ind,:); SX = SX(:,X_inl_ind);
        SY = SY_tmp(Y_inl_ind,:); SY = SY(:,Y_inl_ind);
        
    end
   
    Map_ini = ones(LX,LY)/max(LX,LY);

    
    [Map,fw_output] = FW_minimu_r(Map_ini,M,SX,SY,A,B,opt_gm);
    if opt_gm.r_fixed > 0
        opt_gm.mass_r = sum(Map(:));
    end
    mass1 = min(max(round(sum(Map(:))),1),min(LX,LY));
    
    
    if opt_gm.output > 0
        output.optim_output(cnt).optim_output = fw_output;
        output.Maps(cnt).Maps = Map;
        output.mass_r(cnt) = sum(sum(Map));
    end
 
    if opt_gm.remove_time > 1 && cnt < opt_gm.remove_time
        if cnt <= 1
            th = 0.5;
        else
            th = 0.5;
        end
        
        X_inl_tmp = X_inl_ind;
        Y_inl_tmp = Y_inl_ind;
        
        if ~isfield(opt_gm,'out_re_type')
            [X_inl_ind,Y_inl_ind] = out_re(X_inl_ind,Y_inl_ind,Map,th,mass1);
        else
            [X_inl_ind,Y_inl_ind] = out_re(X_inl_ind,Y_inl_ind,Map,th,mass1,opt_gm.out_re_type);
        end
        
        XX_inl = X(X_inl_ind,:);
        YY_inl = Y(Y_inl_ind,:);
        LX = length(X_inl_ind);
        LY = length(Y_inl_ind);
        
        X_stop = sum(~ismember(X_inl_tmp,X_inl_ind));
        Y_stop = sum(~ismember(Y_inl_tmp,Y_inl_ind));
        stop = (X_stop+Y_stop);
        
        cnt = cnt+1;
        
        if opt_gm.output > 0
            output.inl(cnt).xinl = X_inl_ind;
            output.inl(cnt).yinl = Y_inl_ind;
        end
        
    else 
        stop = 0;
    end
    
    output.X_inl = X_inl_ind;
    output.Y_inl = Y_inl_ind;

end

end