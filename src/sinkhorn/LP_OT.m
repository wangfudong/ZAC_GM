function [Map] = LP_OT(p,q,dis,massrate)
%% solving OT using linear programming
% p,q: marginal constraint
% dis: first order constraint
% massrate: total mass to trasnport
if nargin < 4 || isempty(massrate)
    massrate = 1;
end

LX = length(p);
LY = length(q);

f = dis(:)';

A = zeros(LX+LY,LX*LY);
for i = 1:LY
    A(i,((i-1)*LX+1):i*LX) = 1;
end
for i = 1:LX
    for j = 1:LY
        A(i+LY,(j-1)*LX+i) = 1;
    end
end
b = [q;p];

Aeq = ones(1,LX*LY);
beq = massrate*min(sum(p),sum(q));

lb = zeros(LX*LY,1)';
ub = ones(LX*LY,1)';

tic;
% [X,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub);
[X] = linprog(f,A,b,Aeq,beq,lb,ub);
toc;
Map = reshape(X,[LX,LY]);
% figure,imagesc(Map);title('lp Map');