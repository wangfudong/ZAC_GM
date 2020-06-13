function [Map,output] = ZAC_vgg(X,Y,opt)

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
X_inl_id = (1:LX)';
Y_inl_id = (1:LY)';

if opt.output >0
    output.optim_output = struct([]);
    output.Maps = struct([]);
    output.inl = struct([]);
    output.inl(1).xinl = X_inl_id;
    output.inl(1).yinl = Y_inl_id;
end

mass = opt.mass;
cnt = 1;
stop = 1;

while cnt <= opt.remove_time && stop > 0
   
    DE1 = opt.descx(:,X_inl_id);
    DE2 = opt.descy(:,Y_inl_id);

    M = Dis_sift(DE1',DE2');

    DX = DX_tmp(X_inl_id,:); DX = DX(:,X_inl_id);
    DY = DY_tmp(Y_inl_id,:); DY = DY(:,Y_inl_id);

    SX = SX_tmp(X_inl_id,:); SX = SX(:,X_inl_id);
    SY = SY_tmp(Y_inl_id,:); SY = SY(:,Y_inl_id);
       
    sigx = max(std(DX(DX(:)>0))^2,0.01);
    sigy = max(std(DY(DY(:)>0))^2,0.01);
    sig = 1*[sigx,sigy];
    % sig = [0.5,0.5].^2;
    A = exp(-DX.^2/sig(1)).*(double(DX>0) + eye(size(DX)));
    B = exp(-DY.^2/sig(2)).*(double(DY>0) + eye(size(DY)));
    opt.mass_rate = mass/min(LX,LY);
    
    Map_ini = kLAP_Hun(M,opt.mass);
    %Map_ini = ones(LX,LY)/max(LX,LY);
    
    [Map,fw_output] = FW_minimu(Map_ini,M,SX,SY,A,B,opt);
    
    if opt.output > 0
        output.optim_output(cnt).optim_output = fw_output;
        output.Maps(cnt).Maps = Map;
    end
    
    if opt.remove_time > 1 && cnt < opt.remove_time
        th = 0.5;
        
        X_inl_tmp = X_inl_id;
        Y_inl_tmp = Y_inl_id;
        
        if ~isfield(opt,'out_re_type')
            [X_inl_id,Y_inl_id] = out_re(X_inl_id,Y_inl_id,Map,th,mass);
        else
            [X_inl_id,Y_inl_id] = out_re(X_inl_id,Y_inl_id,Map,th,mass,opt.out_re_type);
        end
        
%         XX_inl = X(X_inl,:);
%         YY_inl = Y(Y_inl,:);
        LX = length(X_inl_id);
        LY = length(Y_inl_id);
        
        X_stop = sum(~ismember(X_inl_tmp,X_inl_id));
        Y_stop = sum(~ismember(Y_inl_tmp,Y_inl_id));
        stop = (X_stop+Y_stop);
        cnt = cnt+1;
        
        if opt.output > 0
            output.inl(cnt).xinl = X_inl_id;
            output.inl(cnt).yinl = Y_inl_id;
        end
        
    else 
        stop = 0;
    end
    
    output.X_inl = X_inl_id;
    output.Y_inl = Y_inl_id;
    
end

end