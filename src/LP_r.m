function [Map] = LP_r(p,q,dis)
%% solving OT using linear programming
% p,q: marginal constraint
% dis: first order constraint
% massrate: total mass to trasnport
% if nargin < 4 || isempty(massrate)
%     massrate = 1;
% end

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

lb = zeros(LX*LY,1)';
ub = ones(LX*LY,1)';

%tic;
% [X,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub);
options.Algorithm = 'interior-point';%'dual-simplex';
options.Display = 'off';
% options.MaxTime = 0.01;
% [X] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
% [X] = linprog(f,A,b,[],[],lb,ub,[],options);
[X] = linprog_fast(f,A,b,[],[],lb,ub,[],options);
%[X] = linprog(f',A,b,[],[],lb',ub',options);%for Gurobi and mosek
%toc;
if isempty(X)
    options.Algorithm = 'dual-simplex';
%     [X] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
    [X] = linprog(f,A,b,[],[],lb,ub,[],options);

end
Map = reshape(X,[LX,LY]);
