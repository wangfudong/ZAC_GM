function [Map_next,output] = FW_minimu(Map_ini,M0,W1,W2,A,B,option)
% Frank-Wolfe method for minimizing the objective function
% M0: node dissimilarity matrix of size m*n
% W1: weight matrix of size m*m
% W2: weight matrix of size n*n
% A: edge attributes matrix of size m*m
% B: edge attributes matrix of size n*n
% option: options of algorithm

max_iter = option.maxiter;
min_iter = option.miniter;

[m,n] = size(M0);

ST = 10;
b = 1/2;
kp = 0;
diff_right = 1;
cnt = 1;
map_tmp = Map_ini;

F_handle = obj_handle(M0,W1,W2,A,B,option);
[FH1, FH2, FH3] = obj_handle_check(M0,W1,W2,A,B,option);

M = obj_gradient(M0,Map_ini,W1,W2,A,B,option);

if option.output > 0
    kp_out = zeros(1,min_iter);
    gra_out = zeros(min_iter,m,n);
    val_out1 = zeros(1,min_iter);
    val_out2 = zeros(1,min_iter);
    val_out3 = zeros(1,min_iter);
    maps_out = zeros(min_iter,m,n);
    maps_out(1,:,:) = Map_ini;
    val_out1(1) = FH1(Map_ini);
    val_out2(1) = FH2(Map_ini);
    val_out3(1) = FH3(Map_ini);
end

if option.active > 1
    active = option.active;
    hist_direc = zeros(m,n,active);
    %weights = ones(1,active)/active;
    weights = exp(-((1:active)-active).^2/active);
    weights = weights/sum(weights);
end

condition1 = 1;
condition2 = 1;
small_tol = 1.0e-5;

while condition1 > 0 || condition2 > 0

    Map = kLAP_Hun(M,option.mass);

    if cnt <= 30% for less iteration
        kp = 0;
    else
        kp = 3;
    end
    
    if option.active > 1% for modified FW method. default: not use
        hist_direc(:,:,active+1) = Map;
        hist_direc(:,:,1) = [];
        dmap = 0;
        for ii = 1:active
            dmap = dmap + weights(ii)*(hist_direc(:,:,ii));
        end
        dmap = dmap - map_tmp;
    else
        dmap = Map - map_tmp;
    end
    
    
    diff_left = -sum(sum(M.*dmap));
    diff_right = 1;
    if abs(diff_left) >= small_tol
        while kp <= ST && abs(diff_right) >= small_tol
            F1 = F_handle(map_tmp);
            F2 = F_handle(map_tmp+b^kp*dmap);
            diff_right = F1 - F2;
            if diff_right >= small_tol
                %break;
                diff_right0 = F1 - F_handle(map_tmp+b^(kp+1)*dmap);
                if diff_right0 > diff_right
                    kp = kp + 1;
                    break;
                else
                    break;
                end
                
            else
                kp = kp + 1;
            end
        end
        
        Map_next = map_tmp + b.^kp*dmap;
        map_tmp = Map_next;
        
    else
        kp = ST + 1;
        Map_next = Map;
    end
    
    if option.output > 0
        kp_out(cnt) = kp;
        gra_out(cnt,:,:) = M;
        val_out1(cnt+1) = FH1(Map_next);
        val_out2(cnt+1) = FH2(Map_next);
        val_out3(cnt+1) = FH3(Map_next);
        maps_out(cnt+1,:,:) = Map_next;
    end

    cnt = cnt +1;
    condition1 = (cnt <= min_iter);
    condition2 = (cnt <= max_iter)*(kp <= ST)*(diff_right >= small_tol);
    
    if  kp <= ST || condition1
        M = obj_gradient(M0,Map_next,W1,W2,A,B,option);
    end
    
end

if option.output > 0
    output.step = kp_out;
    output.value1 = val_out1;
    output.value2 = val_out2;
    output.value3 = val_out3;
    output.gra_out = gra_out;
    output.maps_out = maps_out;
else
    output.step = [];
    output.value = [];
    output.gra_out = [];
    output.maps_out = [];
end

disp(['  Last cg_Iter:',num2str(cnt-1),' Last cg_Err: ',num2str(diff_right) ' k = ' num2str(kp)]);

end





