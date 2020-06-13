function [ X score ] = GAM( W, L12, conf1, conf2 )
% function [ X score ] = GAM( W, L12, conf1, conf2 , par)

W = W - 1.4*eye(size(W));

i = 1;
while 1
    if L12(i,1) == L12(i+1,1)
        order = 'col-wise'; break;
    elseif L12(i,2) == L12(i+1,2)
        order = 'row-wise'; break;
    else
        i = i+1;
    end
    if i == length(L12)-1
        order = 'none'; break;
    end
end
if strcmp(order, 'col-wise')
    E12 = accumarray(L12, 1)';
elseif strcmp(order, 'row-wise')
    E12 = accumarray(L12, 1);
end
[n1,n2]=size(E12);

% b0: sensitive parameter to tweak
coeff=1;
b0=coeff*max(n1,n2);

% bMax=1e3;
bMax=500;  % ECCV 2012
% bMax = par;
% bMax = 1700;

%tolB=1e-10;%1e-3;
tolB=1e-8;%1e-3;
% tolC=1e-3;%1e-3;
tolC=1e-3;%1e-3;

bStep = 1.075; % ECCV 2012
% bStep = par;
% bStep = 1.05;

W = sparse(W);
%[X2,nbMatVec] = gradAssign(W, E12, b0, 1.075,bMax,tolB, tolC, target,X);
[X2,nbMatVec] = gradAssign(W, E12, b0, bStep,bMax,tolB, tolC, conf1, conf2 ); % ECCV 2012
% [X2,nbMatVec] = gradAssign(W, E12, 100, 1,bMax,tolB, tolC, conf1, conf2 ); % PR project

%% post processing for scoring 
X = X2(find(E12));

