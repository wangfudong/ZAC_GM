function [KP, KQ] = conKnlGphPQU_tmp(gphs, parKnl)
% Compute node and feature affinity matrix for graph matching.
%
% Remarks
%   The edge is undirected.
%   To deal with directed edge, please use the function "conKnlGphPQD.m".
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%   parKnl  -  parameter
%     alg   -  method of computing affinity, {'toy'} | 'cmum' | 'pas'
%              'toy':  toy data
%              'cmum': CMU motion data
%              'pas':  Pascal data
%
% Output
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% function parameter
alg = ps(parKnl, 'alg', 'toy');

% dimension
gph1 = gphs{1};
gph2 = gphs{2};
[n1, m1] = size(gph1.G);
[n2, m2] = size(gph2.G);
m1 = m1 * 2;
m2 = m2 * 2;
prIn('conKnlGphPQU', 'alg %s, n1 %d, n2 %d, m1 %d, m2 %d', alg, n1, n2, m1, m2);

% for toy data
if strcmp(alg, 'toy')
    KP = zeros(n1, n2);
    DQ = conDst(gph1.XQ, gph2.XQ);
    KQ = exp(-DQ / .15);
    
    % for CMU Motion data
elseif strcmp(alg, 'cmum')
    KP = zeros(n1, n2);
    
%     DP = gphs{1}.sc;
%     KP = exp(-DP/1);
    DQ = conDst(gph1.dsts, gph2.dsts);
    KQ = exp(-DQ / 2500);
    
    % for point registration
elseif strcmp(alg, 'pr')
    DP = gphs{1}.unary;
    KP = exp(-DP);
    KP = zeros(n1,n2);
    DQ = conDst(gph1.dsts, gph2.dsts);
    KQ = exp(-DQ / 2500);
    
    % for Pascal data
elseif strcmp(alg, 'pas')
    
    DP = conDst(gphs{1}.Pt, gphs{2}.Pt);
    KP = exp(-real(sqrt(DP)));
    
%     DP = gphs{1}.sc;
%     KP = exp(-DP/0.1);
    
    % distance
    Dst1 = repmat(gph1.dsts', 1, m2);
    Dst2 = repmat(gph2.dsts, m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);
    
    % angle
    Ang1 = repmat(gph1.angs', 1, m2);
    Ang2 = repmat(gph2.angs, m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    
elseif strcmp(alg, 'face205')
%     [~,~,hks1,hks2] = get_graph_face205(gphs{1}.orders);
    hks1 = gphs{1}.hks1;
    hks2 = gphs{2}.hks2;
    DP = measure_hks(hks1,hks2,'E');
    KP = exp(-real(DP));
    
    gind1 = gph1.Eg;
    gind2 = gph2.Eg;
    dsts1 = zeros(1,length(gind1));
    dsts2 = zeros(1,length(gind2));
    for i = 1:length(gind1)
        dsts1(i) = gphs{1}.DX(gind1(1,i),gind1(2,i));
    end
    for ii = 1:length(gind2)
        dsts2(ii) = gphs{2}.DY(gind2(1,ii),gind2(2,ii));
    end
    
    DQ = conDst(dsts1,dsts2);
    KQ = exp(-DQ / 500);
    
%     % distance
%     Dst1 = repmat(gph1.dsts', 1, m2);
%     Dst2 = repmat(gph2.dsts, m1, 1);
%     Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);
%     
%     % angle
%     Ang1 = repmat(gph1.angs', 1, m2);
%     Ang2 = repmat(gph2.angs, m1, 1);
%     Ang = abs(Ang1 - Ang2);
%     
%     % combine distance and angle
%     KQ = exp(-(Dst + Ang) / 2);
    
else
    error('unknown algorithm: %s', alg);
end

% normalize
% KQ = knlEgNor(KQ, parKnl);

% for symmetric edges
KQ = KQ(1 : round(m1 / 2), 1 : round(m2 / 2));

prOut;
