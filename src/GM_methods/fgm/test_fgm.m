%% test fgmdemos
% clear variables; 
% prSet(1);
%%
dataset = ['house';'hotel'];
sizeimage = [576,512];
data = 2;g1 = 1;
numb1 = num2str(g1,'%03d');
file1 = [dataset(data,:) numb1];
hl1 = load(file1);
hi1 = imread([dataset(data,:) '.seq' num2str(g1-1) '.png']);
LX = 20;
LY = 30;

  order = randperm(30);
  order = order(1:LX);
% % % %order = 1:1:LX;

g2 = 50;
numb2 = num2str(g2,'%03d');file2 = [dataset(data,:) numb2];
hl2 = load(file2);
hi2 = imread([dataset(data,:) '.seq' num2str(g2-1) '.png']);

XX = hl1(order,:);
YY = hl2;

Fs{1} = hi1;
Fs{2} = hi2;
Pts{1} = XX';
Pts{2} = YY';
%% algorithm parameter
gt = zeros(LX,LY);
for i = 1:LX
    gt(i,order(i)) = 1;
end
asgT.alg = 'truth';
asgT.X = gt;

parKnl = st('alg', 'pas'); % type of affinity: only edge distance
[pars, algs] = gmPar(2);
%% src
%wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
%asgT = wsSrc.asgT;

%% feature
parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
gphs = newGphUs(Pts, parG);

%% affinity
[KP, KQ] = conKnlGphPQU(gphs, parKnl);%%parKnl 可以用 pas 表示Pascal 数据集
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

%% undirected graph -> directed graph (for FGM-D)
gphDs = gphU2Ds(gphs);
KQD = [KQ, KQ; KQ, KQ];

%% GA
asgGa = gm(K, Ct, asgT, pars{1}{:});

%% PM
asgPm = pm(K, KQ, gphs, asgT);

%% SM
asgSm = gm(K, Ct, asgT, pars{3}{:});

%% SMAC
asgSmac = gm(K, Ct, asgT, pars{4}{:});

%% IPFP-U
asgIpfpU = gm(K, Ct, asgT, pars{5}{:});

%% IPFP-S
asgIpfpS = gm(K, Ct, asgT, pars{6}{:});

%% RRWM
asgRrwm = gm(K, Ct, asgT, pars{7}{:});

%% FGM-U
asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:});

%% FGM-D
asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, pars{9}{:});

%% print information
fprintf('GA    : acc %.2f, obj %.2f\n', asgGa.acc, asgGa.obj);
fprintf('PM    : acc %.2f, obj %.2f\n', asgPm.acc, asgPm.obj);
fprintf('SM    : acc %.2f, obj %.2f\n', asgSm.acc, asgSm.obj);
fprintf('SMAC  : acc %.2f, obj %.2f\n', asgSmac.acc, asgSmac.obj);
fprintf('IPFP-U: acc %.2f, obj %.2f\n', asgIpfpU.acc, asgIpfpU.obj);
fprintf('IPFP-S: acc %.2f, obj %.2f\n', asgIpfpS.acc, asgIpfpS.obj);
fprintf('RRWM  : acc %.2f, obj %.2f\n', asgRrwm.acc, asgRrwm.obj);
fprintf('FGM-U : acc %.2f, obj %.2f\n', asgFgmU.acc, asgFgmU.obj);
fprintf('FGM-D : acc %.2f, obj %.2f\n', asgFgmD.acc, asgFgmD.obj);
%% show correspondence result
rows = 1; cols = 1;
Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b', 'g'});
shAsgImg(Fs, gphs, asgFgmD, asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
title('result of FGM-D');

%% show affinity
rows = 1; cols = 3;
Ax = iniAx(2, rows, cols, [200 * rows, 200 * cols]);
shAsgK(K, KP, KQ, Ax);

%% show correpsondence matrix
asgs = {asgT, asgGa, asgPm, asgSm, asgSmac, asgIpfpU, asgIpfpS, asgRrwm, asgFgmU, asgFgmD};
rows = 2; cols = 5;
Ax = iniAx(3, rows, cols, [250 * rows, 250 * cols]);
shAsgX(asgs, Ax, ['Truth', algs, 'FGM-A']);
