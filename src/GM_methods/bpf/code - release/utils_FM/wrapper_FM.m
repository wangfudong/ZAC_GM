%% Makes current problem into Graph Matching form
% Note: 'cdata' should contain all variables that GM solver needs
%       ex) RRWM solver needs: affineMatrix, group1, group2
function [accuracy score time X Xraw] = wrapper_FM(method, cdata)
% Make function evaluation script
str = ['feval(@' func2str(method.fhandle)];
for j = 1:length(method.variable), str = [str ',cdata.' method.variable{j} ]; end
if ~isempty(method.param), for i = 1:length(method.param), str = [str, ',method.param{' num2str(i) '}']; end; end
str = [str, ')']; 
% Function evaluation & Excution time Check
tic; Xraw = eval(str); time = toc;
X = greedyMapping(Xraw, cdata.group1, cdata.group2);
% Matching Score
score = X'*cdata.affinityMatrix*X; % objective score function

% extrapolate the solution for flexible evaluation
X = extrapolateMatchIndicator(cdata.view, cell2mat({cdata.matchInfo.match}'),X,cdata.extrapolate_thres)';
if length(cdata.GTbool) ~= length(cdata.affinityMatrix)
    accuracy = NaN; % Exception for no GT information
else
    accuracy = (X(:)'*cdata.GTbool(:))/sum(cdata.GTbool);
end