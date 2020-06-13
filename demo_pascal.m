%% demo of pascal dataset
dataset = ['Carss';'Motor'];
dataL = [30,20];% number of images in the dataset

datasets = 3;%3=Carss,4=Motors
inds = 9;%index of graph pairs 

[Pts,Ims,nF] = fileload(datasets,inds);
X_all = Pts{1,1};
Y_all = Pts{1,2};
I1 = Ims{1,1};
I2 = Ims{1,2};

max_size1 = max(size(I1));max_size2 = max(size(I2));
max_size = max(max_size1,max_size2);

[X_inl,X_out] = pascal_pick_outlier(X_all,1:nF,[],20);
[Y_inl,Y_out] = pascal_pick_outlier(Y_all,1:nF,[],20);

gap = 10;
pad_size = [400,400,3];
[Xnew,Ynew,im12] = pascal_image_padding(I1,I2,[X_inl;X_out],[Y_inl;Y_out],gap,pad_size);

shifted = 0;
x_shift = pad_size(2) + gap;
hh=figure;set(hh,'position',[500,300,700,430]);
imshow(im12);hold on;
set(gca,'position',[0.01,0.01,0.98,0.99]);
plot(Xnew(1:nF,1),Xnew(1:nF,2)-shifted,'r.','markersize',20);hold on;
%text(Xnew(1:nF,1),Xnew(1:nF,2)-shift1,num2str([1:nF]'),'color','g','fontsize',13);hold on;
plot(Xnew(nF+1:end,1),Xnew(nF+1:end,2)-shifted,'y*','markersize',10,'linewidth',1.5);hold on;
plot(Ynew(1:nF,1)+x_shift,Ynew(1:nF,2)-shifted,'r.','markersize',20);hold on;
%text(Ynew(1:nF,1)+x_shift,Ynew(1:nF,2)-shift1,num2str([1:nF]'),'color','g','fontsize',13);hold on;
plot(Ynew(nF+1:end,1)+x_shift,Ynew(nF+1:end,2)-shifted,'y+','markersize',10,'linewidth',2.5);hold on;
drawnow;

order = 1:nF;
%order = randperm(nF);
Xtmp = Xnew(1:nF,:);
XX = [Xtmp(order,:);Xnew(nF+1:end,:)]/max_size;
YY = Ynew/max_size;
%% ZAC
GT = [1:length(order);order]';
opt = set_option('w/o',nF,GT);

tic;
[Map,output] = ZAC(XX,YY,opt);
toc;

% [re,pre,xinlrate,yinlrate] = acc_compute(Map,GT,output.X_inl,output.Y_inl,nF);
plotMatch_pascal(asgHun(Map),GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],output.X_inl,output.Y_inl,im12,nF,x_shift);

%% ZACR
GT = [1:length(order);order]';
opt = set_option('w',nF,GT);

tic;
[Map,output] = ZAC_r(XX,YY,opt);
toc;

% [re,pre,xinlrate,yinlrate] = acc_compute(Map,GT,output.X_inl,output.Y_inl,nF);
plotMatch_pascal(asgHun(Map),GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],output.X_inl,output.Y_inl,im12,nF,x_shift);


%% the other GM methods
Pts{1,1} = (XX*max_size)';
Pts{1,2} = (YY*max_size)';

parKnl = st('alg', 'pas'); % type of affinity: only edge distance
parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
gphs = newGphUs(Pts, parG);

[~, KQ] = conKnlGphPQU(gphs, parKnl);
KP =  M_shape(Pts{1,1}',Pts{1,2}',1/8,2,0);
KP = exp(-KP/1);
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

[ind, ind_m] = find(Ct);
[m,n] = size(KP);
group1 = zeros(size(ind, 1), m);
group2 = zeros(size(ind, 1), n);
for i = 1:size(ind, 1)
    group1(i, ind(i)) = 1;
    group2(i, ind_m(i)) = 1;
end
group1 = logical(group1);
group2 = logical(group2);
gt = zeros(length(XX(:,1)),length(YY(:,1)));
for i = 1:nF
    gt(i,order(i)) = 1;
end
asgT.X = gt;
%% GAGM
t1 = clock;
E12 = ones(size(group1, 2), size(group2,2));
[L1, L2] = find(E12);
L12 = [L1, L2];
GAM_soft = GAM(K, L12, group1, group2);
GAM_bin = greedyMapping(GAM_soft, group1, group2);
GAM_bin = reshape(GAM_bin, m, n);
t2 = clock;
asgGagm.tm = etime(t2,t1);
asgGagm.acc = matchAsg(GAM_bin, asgT);
plotMatch_pascal(GAM_bin,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)

GAM_topK = Map_by_TopK(reshape(GAM_soft, m, n),GAM_bin,nF);
asgGagm.acc = matchAsg(GAM_topK, asgT);
plotMatch_pascal(GAM_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)

%% IPFP-S
t1 = clock;
IPFPS_soft = IPFP_S(K, group1, group2);
IPFPS_bin = greedyMapping(IPFPS_soft, group1, group2);
IPFPS_bin = reshape(IPFPS_bin, m, n);
t2 = clock;
asgIpfpS.tm = etime(t2,t1);
asgIpfpS.acc = matchAsg(IPFPS_bin, asgT);
plotMatch_pascal(IPFPS_bin,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);

IPFPS_topK = Map_by_TopK(reshape(IPFPS_soft, m, n),IPFPS_bin,nF);
asgGagm.acc = matchAsg(IPFPS_topK, asgT);
plotMatch_pascal(IPFPS_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);
%% PSM
t1 = clock;
PSM_soft = PSM(K, group1, group2);
PSM_bin = greedyMapping(PSM_soft, group1, group2);
PSM_bin = reshape(PSM_bin, m, n);
t2 = clock;
asgPsm.tm = etime(t2,t1);
asgPsm.acc = matchAsg(PSM_bin, asgT);
plotMatch_pascal(PSM_bin,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);

PSM_topK = Map_by_TopK(reshape(PSM_soft, m, n),PSM_bin,nF);
asgPsm.acc = matchAsg(PSM_topK, asgT);
plotMatch_pascal(PSM_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);
%% RRWM
t1 = clock;
RRWM_soft = RRWM(K, group1, group2);
RRWM_bin = greedyMapping(RRWM_soft, group1, group2);
RRWM_bin = reshape(RRWM_bin, m, n);
t2 = clock;
asgRrwm.tm = etime(t2,t1);
asgRrwm.acc = matchAsg(RRWM_bin, asgT);

plotMatch_pascal(RRWM_bin,GT,[XX(:,1)*max_size,XX(:,2)*max_size],[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)
% [re,pre,xinlrate,yinlrate] = acc_compute(RRWM_bin,GT,(1:m)',(1:n)',nF);

RRWM_topK = Map_by_TopK(reshape(RRWM_soft, m, n),RRWM_bin,nF);
% [re,pre,xinlrate,yinlrate] = acc_compute(RRWM_topK,GT,(1:m)',(1:n)',nF);
plotMatch_pascal(RRWM_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)
%% BPF_G
t1 = clock;
parGnccp = st('nItMa', 100, 'deta', 0.001, 'nHist', 5, 'rho', 2, 'theta', 0.01);
X = BPF_G(K, group1, group2, parGnccp);
X = greedyMapping(X, group1, group2);
X = reshape(X, m, n);
t2 = clock;
asgBpfG.tm = etime(t2,t1);
asgBpfG.acc = matchAsg(X, asgT);

plotMatch_pascal(X,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)

%% FGM-D
t1 = clock;
parFgmA = st('nItMa', 100, 'nAlp', 101, 'deb', 'n', 'ip', 'n', 'lamQ', .5);
gphDs = gphU2Ds(gphs);
KQD = [KQ, KQ; KQ, KQ];
asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, parFgmA);
t2 = clock;
asgFgmD.tm = etime(t2,t1);

plotMatch_pascal(asgFgmD.X,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)
%% MPM
t1 = clock;
asgMPM = MPM_GM( XX*max_size, YY*max_size,(1:m)',K,2500);
t2 = clock;
asgMPM.tm = etime(t2,t1);
asgMPM.acc = matchAsg(asgHun(asgMPM.X), asgT);
plotMatch_pascal(asgHun(asgMPM.X),GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);

MPM_topK = Map_by_TopK(reshape(asgMPM.X_con, m, n),asgHun(asgMPM.X),nF);
plotMatch_pascal(MPM_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift);

%% FRGM without outlier-removal
option.M_exist = 10;
option.alpha2 = 1;
option.maxiter = 100;
option.active = 1;
option.func = 'inner';
option.q_norm = 2;
option.type = 'inner';
option.full = 1;
option.rota = 0;

[M,~,~] = M_shape(XX,YY,1/8,2,option.rota);
DX = M_points(XX,XX);
DY = M_points(YY,YY);
max_fact = max([max(DX(:)),max(DY(:))]);
DXX = DX/max_fact;DYY = DY/max_fact;
SX = 1./(DXX + 1*(DXX==0)).*(DXX > 0);
SX = SX/max(SX(:));

Map_ini = asgHun(-M);
t1 = clock;
sig_factor = 1;
sig = sig_factor*[std(DXX(:)).^2,std(DYY(:)).^2];
weight = [0.5,0.5];
[asgFRGM] = FRGM_Gen(Map_ini,M,DXX,DYY,SX,option,asgT,sig,weight);
t2 = clock;
plotMatch_pascal(asgHun(asgFRGM.map2),GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)

Map_topK = Map_by_TopK(asgFRGM.map2,asgHun(asgFRGM.map2),nF);
plotMatch_pascal(Map_topK,GT,XX*max_size,[Ynew(:,1),Ynew(:,2)-shifted],(1:m)',(1:n)',im12,nF,x_shift)


