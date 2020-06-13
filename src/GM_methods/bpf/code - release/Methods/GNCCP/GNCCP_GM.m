function [distance, X] = GNCCP_GM(Ct, Ag, Ah, para)
%graph matching by GNCCP
%min. tr(Ct'X) + tr(Ag-XAhX')'(Ag-XAhX')
%Ag, Ah are adjacency matrix, Ng <= Nh
%
% Reference
%   [1] Zhi-Yong Liu and Hong Qiao, "GNCCP - Graduated NonCovexity and Concavity
%       Procedre", IEEE Trans. PAMI, DOI: 10.1109/TPAMI.2013.223, 2014
%   [2] Zhi-Yong Liu, H. Qiao, X. Yang and Steven C.H. Hoi, "Graph Maching by
%       Simplified Convex-Concave Relaxation Procedure", IJCV, DOI: 10.1007/s11263-014-0707-7, 2014
%
% History
%   written by Prof. Zhi-Yong Liu, at Mar. 2014
%   modified by Tao Wang (twang@bjtu.edu.cn),   01-07-2015

nh = size(Ah, 1);
ng = size(Ag, 1);
if ng > nh
    error('input error');
end

P0 = ones(ng, nh)/nh;

%parameter
eta = 0.001; %for gradient value
P = P0;

nItMa = 0;
step = 0.002;
if nargin > 3
    nItMa = para.nItMa;
    step = para.step;
end


gamma = 1;
%t=0
while gamma > -1
    ite = 0;
    while 1
        ite = ite + 1;
        if nItMa > 0 && ite > nItMa
            break;
        end
        
        if gamma > 0
            g = gamma * P + (1-gamma)*(0.5 * Ct + P*(Ah'*P'*P*Ah + Ah*P'*P*Ah') - (Ag*P*Ah' + Ag'*P*Ah));
        else
            g = gamma * P + (1+gamma)*(0.5 * Ct + P*(Ah'*P'*P*Ah + Ah*P'*P*Ah') - (Ag*P*Ah' + Ag'*P*Ah));
        end
%         g = g + Ct;

        if trace(g * g') < eta
            break;
        end
        

        %subproblem 
        G = zeros(nh, nh);
        G(1:ng, :) = g;
        maxv = max(max(g));
        if nh > ng
            G(ng+1:nh, :) = maxv*ones(nh-ng,nh);
        end

        X = KM(-G);
        X = X(1:ng,:);
       
        errorterm = sum(sum(g.*(X-P)));
        if(errorterm > -eta) 
            break; 
        end

        %line search
        %divided method 1
        Pright = X;
        Pleft = P;
      %  FX = fl(Ct, Ah, Ag, X, gamma);
      %  FP = fl(Ct, Ah, Ag, P, gamma);
        deltaX = (Pright - Pleft);
        for i = 1:10
            Pnew = Pleft + 0.5 * deltaX;
            Pnew2 = Pleft + (0.5 + eta) * deltaX;
            F0 = fl(Ct, Ah, Ag, Pnew, gamma);
            F1 = fl(Ct, Ah, Ag, Pnew2, gamma);

            if F0 < F1
               Pright = Pnew;
            else
               Pleft = Pnew;
            end
            deltaX = Pright - Pleft;
        end
        P = Pnew;
        
        Fnew = fl(Ct, Ah, Ag, P, gamma);
        if gamma > 0
            g = gamma * P + (1-gamma)*(0.5 * Ct + P*(Ah'*P'*P*Ah + Ah*P'*P*Ah') - (Ag*P*Ah' + Ag'*P*Ah));
        else
            g = gamma * P + (1+gamma)*(0.5 * Ct + P*(Ah'*P'*P*Ah + Ah*P'*P*Ah') - (Ag*P*Ah' + Ag'*P*Ah));
        end
%         g = g + Ct;
        
        tmp = sum(sum(g.*(P - X)));
        if tmp < eta  * abs(Fnew - tmp);
            break;
        end

    end %end of frank-wolfe
    if nh - trace(P'*P) < 0.1^3
        break;
    end
    gamma = gamma - step;
    if mod(round(gamma * 1000), 100) == 0
        disp(['GNCCP_GM: gamma = ', num2str(gamma)]);
    end
end %end of gamma

PG = zeros(nh, nh);
PG(1:ng, :) = P;
minv = min(min(P));
if nh > ng
    PG(ng+1:nh, :) = minv*ones(nh-ng,nh);
end
X = KM(PG);
X = X(1:ng,:);

distance = trace(Ct'*X) + sum(sum((Ag - X * Ah * X') .* (Ag - X * Ah * X')));
end

function f = fl(Ct, Ah, Ag, P, gamma)
    if gamma > 0
        f = gamma * sum(sum((P.*P))) + (1-gamma) * (sum(sum(Ct .* P)) + sum(sum((Ag - P * Ah * P').*(Ag - P * Ah * P'))) );
    else
        f = gamma * sum(sum((P.*P))) + (1+gamma) * (sum(sum(Ct .* P)) + sum(sum((Ag - P * Ah * P').*(Ag - P * Ah * P'))) );
    end
   
end
