%% non-rigid deformation example
data_name = 5;%1=beijing;2=whale;3=chinese;4=cpd_fish;5=fish_2d;
X0 = load_testdata(data_name);
LX0 = length(X0(:,1));
[X0,X_bound] = normalize_point(X0,1);

sample_rate = 0.7;
samples = randperm(LX0);
LX_sample = round(sample_rate*LX0);
samples = samples(1:LX_sample);
X = X0(samples,:);

sigma = 2;
GX = nonrigid_ker(X,sigma,'rbf');
theta = 0.1*rand(1)*pi;
R = [cos(theta),sin(theta);
    -sin(theta),cos(theta)];
W = 0.5*randn(length(X(:,1)),2);
t = 1*rand(1,2) + 0.5;
%t=[0,0];
scale = 1;

opt_non.noise = 1;
opt_non.noise_sigma = 0.0;
opt_non.noise_type = 'uniform';%'guassian'
opt_non.outlier = 50;
opt_non.out_sigma = 0.5;
opt_non.outlier_type = 'gaussian';%'uniform1''guassian''guassian1'

Y = scale*nonrigid_kernel_trans(X,W,GX,R,opt_non);
LY = length(Y(:,1));
Y = Y - repmat(mean(Y),LY,1) + repmat(t,LY,1);
opt_non.outlier = 50;
X = rigid_affine_transform2D(X,eye(2,2),[0,0],1,opt_non);

figure,plot(X(1:LX_sample,1),X(1:LX_sample,2),'r.',X(LX_sample+1:end,1),X(LX_sample+1:end,2),'g.','markersize',20);hold on;
plot(Y(:,1),Y(:,2),'g.',Y(1:LX_sample,1),Y(1:LX_sample,2),'b.','markersize',20);axis equal


order = randperm(LX_sample);
%order = 1:LX_sample;
X_inl = X(1:LX_sample,:);
X(1:LX_sample,:) = X_inl(order,:);

gt = zeros(LX_sample,LX_sample);
for i = 1:LX_sample
    gt(i,order(i)) = 1;
end
asgT.X = gt;
LX = length(X(:,1));
gt_all = zeros(LX,LY);
for i = 1:LX_sample
    gt_all(i,order(i)) = 1;
end
%% ZAR and ZACR
nF = LX_sample;
% match_rate = 0.5;
%nF = round(match_rate*LX)
% nF = nF_by_kmeans(X,Y,0,0)

GT = [1:length(order);order]';

opt_gm = set_option('w/o',nF,GT);
% opt_gm = set_option('w',nF,GT);%for ZACR
opt_gm.unary = 10;
opt_gm.remove_time = 5;
opt_gm.output = 0;


opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.reg_maxiter = 30;
opt_reg.rota = 0;
opt_reg.reg_normalize = 1;
opt_reg.nonrigid_sig = 2;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];
opt_reg.lambda = 0.1;
opt_reg.type = 'nonrigid';

tic;
para = ZAC_nonrigid_2d(X,Y,opt_reg,opt_gm);
% para = ZAC_nonrigid_2d_r(X,Y,opt_reg,opt_gm);%for ZACR

toc;
% X_trans = para.sy0*para.X_new + repmat(para.uy,LX,1);
% X_tmp = para.sx0*para.X_tmp + repmat(para.ux,LX,1);
% Y_tmp = para.sy0*para.Y_tmp + repmat(para.uy,LY,1);

figure,plot(1:length(para.reg_Err),para.reg_Err,'r.-');hold on;
plot(1:length(para.reg_Diff),para.reg_Diff,'b.-');
YLIMmax = max(max(para.reg_Err),max(para.reg_Diff));
set(gca,'YLim',[0,YLIMmax+0.05]);

%% majiayi tip
opt.outliers = 1-LX_sample/LX;
tic;
[Transform1, C] = GLS_ma(X,Y,100,1,1);
toc;
reg_err_gls = measurement_XY(Transform1.Y,Y,Y,opt_gm.GT);
%% CPD
opt.method='nonrigid';
opt.corresp = 1;      % compute correspondence vector at the end of registration (not being estimated by default)
opt.outliers = 1-LX_sample/LX;
opt.viz = 1;
opt.inl_num = LX_sample;
[Transform,C_cpd]=cpd_register(Y,X,opt);

% figure,cpd_plot_iter(XT, X); title('Before');
% figure,cpd_plot_iter(Transform.Y,XT);  title('After');

figure,plot(Transform.Y(:,1),Transform.Y(:,2),'bo',Y(:,1),Y(:,2),'r.');

reg_err_cpd = measurement_XY(Transform.Y,Y,Y,opt_gm.GT)


%% GAGM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];
opt_reg.nonrigid_sig = 2;
opt_reg.lambda = 1;

para_gagm = GAGM_nonrigid_2d(X,Y,opt_reg,opt_gm);

%% IFFPS
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];
opt_reg.nonrigid_sig = 2;
opt_reg.lambda = 1;

para_ipfps = IPFPS_nonrigid_2d(X,Y,opt_reg,opt_gm);

%% PSM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];
opt_reg.nonrigid_sig = 2;
opt_reg.lambda = 1;

para_psm = PSM_nonrigid_2d(X,Y,opt_reg,opt_gm);
%% RRWM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];
opt_reg.nonrigid_sig = 2;
opt_reg.lambda = 1;

para_rrwm = RRWM_nonrigid_2d(X,Y,opt_reg,opt_gm);
%% FGMD
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_fgmd = FGMD_nonrigid_2d(X,Y,opt_reg,opt_gm);

%% MPM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_mpm = MPM_nonrigid_2d(X,Y,opt_reg,opt_gm);

%% FRGM
opt_gm = set_opt_FRGM(10,1,1);
opt_gm.GT = GT;
opt_gm.topk = LX_sample;
opt_gm.tol_cg = 1.0e-6; 
opt_gm.Mexist = 1;
opt_gm.maxiter_nonvex = 100;
opt_gm.maxiter_convex = 100;
opt_gm.lam_nonvex = 1;
opt_gm.lam_convex = 1;
opt_gm.lam_sparse = 0;
opt_gm.geofunc = '1.21';
opt_gm.unary = 1;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 0;
opt_reg.reg_maxiter = 30;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_frgm = FRGM_nonrigid_2d_fast(X,Y,opt_reg,opt_gm);









