function X = IPFP(K, x0, group1, group2)
% This function tries to maximize the matching score x' M x 
% where x obeys discrete one-to-one matching constraints 
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
%
% Math
%   This algorithm is to obtain the optimal X for the following problem
%     max_x   x' * K * x
%     s.t.    A * x <= 1
%
% Input
%   K:      [m * m] affinity matrix of m candicate matches
%   x0:     initial assignment vector
%   group1: conflicting match groups in domain 1 (size(K,1) x nGroup1)
%   group2: conflicting match groups in domain 2 (size(K,1) x nGroup2)
%                 
%       e.g. find(group1(:,3)) represents the third goup of matches  
%                               sharing the same feature in domain1   
%
% Output
%   X        -  assignment vector X
%
% History
%   create   -  Marius Leordeanu (leordeanu@gmail.com), 02-25-2011
%   modify   -  Tao Wang (twang@bjtu.edu.cn), 01-03-2015

% function parameter
nItMa = 50;
isDeb = 0;

M = K;

% initial
new_sol = x0;
best_sol = x0;
best_score = 0; 
nSteps = 0;
scores(1) = new_sol' * M * new_sol;
scores2(1) = new_sol' * M * new_sol;
discreteRate = 0;
while nSteps <= nItMa
   nSteps = nSteps + 1;

   old_sol = new_sol;
   xx = M * old_sol;
   
   % projection to discrete domain
   x2 = greedyMapping(xx, group1, group2);

   % step size
   D = (x2 - old_sol)' * M * (x2 - old_sol);
   if D >= 0
       new_sol = x2;
       stepSize_t(nSteps) = 1;
       stepSize_norm(nSteps) = norm(x2 - old_sol);
       discreteRate = discreteRate + 1;
   else
       C = old_sol' * M * (x2 - old_sol);
       r = min([1, -C / D]);
       if r < 0.01
           r = 0;
       elseif r == 1
          discreteRate = discreteRate + 1;
       end
       new_sol = old_sol + r * (x2 - old_sol);
       stepSize_t(nSteps) = r;
       stepSize_norm(nSteps) = norm(x2 - old_sol);
   end

   scores(nSteps + 1) = new_sol' * M * new_sol;
   scores2(nSteps + 1) = new_sol' * M * old_sol;
   dX(nSteps) = sum(abs(new_sol - best_sol));
   curr_score = x2' * M * x2;
   
   if curr_score > best_score
        best_score = curr_score;
        best_sol = x2;
   end
   
   % stop condition
   if norm(new_sol - old_sol) == 0
       break;
   end
end

discreteRate = discreteRate / nSteps;
sol = best_sol;
stats.dX = dX;
stats.scores = scores;
stats.scores2 = scores2;
stats.best_score = best_score;
stats.discreteRate = discreteRate;
stats.stepSize_t = stepSize_t;
stats.stepSize_norm = stepSize_norm;

X = sol;
end

