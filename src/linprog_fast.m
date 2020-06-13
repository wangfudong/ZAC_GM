function [x, fval, exitflag, output, lambda] = linprog_fast(f,A,B,Aeq,Beq,lb,ub,algoptions)

% Alg = options.Algorithm;
% Dis = options.Display;
        
%algoptions = optimoptions('linprog', 'Algorithm', Alg);

% optFcn = optim.options.Linprog;
% algoptions = optFcn('Algorithm', Alg);

% algoptions = optim.options.Linprog;
% algoptions.Display = Dis;
% algoptions.Algorithm = Alg;


nvars = length(f);
if isempty(A), A=zeros(0,nvars); end
if isempty(B), B=zeros(0,1); end
if isempty(Aeq), Aeq=zeros(0,nvars); end
if isempty(Beq), Beq=zeros(0,1); end

problem.f = f;
problem.Aineq = A;
problem.bineq = B;
problem.Aeq = Aeq;
problem.beq = Beq;
problem.lb = lb;
problem.ub = ub;
problem.options = algoptions;
problem.solver = 'linprog';

%algorithm = createAlgorithm(algoptions);
algorithm = optim.algorithm.LinprogInteriorPoint(algoptions);

[x, fval, exitflag, output, lambda] = run(algorithm, problem);

end




