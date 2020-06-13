% demo for point registration :rigid, some fast algorithms
data_name = 5;%1=beijing;2=whale;3=chinese;4=cpd_fish;5=fish_2d;
X0 = load_testdata(data_name);
LX0 = length(X0(:,1));
[X0,X_bound] = normalize_point(X0,1);

samples = randperm(LX0);
LX_sample = round(0.8*LX0);
samples = samples(1:LX_sample);
X = X0(samples,:);

opt_similar.noise = 1;
opt_similar.noise_sigma = 0.01;
opt_similar.noise_type = 'uniform';%or 'gaussian'
opt_similar.outlier = 50;
opt_similar.out_sigma = 0.5;
opt_similar.outlier_type = 'gaussian';%or 'uniform1''gaussian''gaussian1'
    
theta = rand(1)*sign(rand-0.5)*pi;
%theta = -0.1*pi;
R = [cos(theta),sin(theta);
    -sin(theta),cos(theta)];
t = 0.5*X_bound*rand(1,2) + 0.1;
s = 0.1;

Y = rigid_affine_transform2D(X,R,t,s,opt_similar);% add deformation, noise and outliers
X = rigid_affine_transform2D(X,eye(2,2),[0,0],1,opt_similar);

figure();
subplot(1,2,1),plot(X(1:LX_sample,1),X(1:LX_sample,2),'r.',X(LX_sample+1:end,1),X(LX_sample+1:end,2),'g.','markersize',20);
subplot(1,2,2),plot(Y(:,1),Y(:,2),'g.',Y(1:LX_sample,1),Y(1:LX_sample,2),'b.','markersize',20);

LY = length(Y(:,1));
order = randperm(LX_sample);
%order = 1:LX;
X_inl = X(1:LX_sample,:);
X(1:LX_sample,:) = X_inl(order,:);

gt = zeros(LX_sample,LX_sample);
for i = 1:LX_sample
    gt(i,order(i)) = 1;
end
asgT.X = gt;
LX = length(X(:,1));
%% ZAC and ZACR
% nF = LX_sample;
nF = round(0.5*LX);
GT = [1:length(order);order]';

opt_gm = set_option('w/o',nF,GT);
% opt_gm = set_option('w/o',nF,GT);% for ZACR
opt_gm.unary = 10;
opt_gm.remove_time = 5;
opt_gm.output = 0;

opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.reg_maxiter = 20;
opt_reg.rota = 1;
opt_reg.reg_normalize = 1;
opt_reg.unary_terms = [zeros(1,3),ones(1,20)];
opt_reg.type = 'rigid';

t11 = clock;
para_our = ZAC_similar_2d(X,Y,opt_reg,opt_gm);
% para_our = ZAC_similar_2d_r(X,Y,opt_reg,opt_gm);% for ZACR

t22 = clock;
time12 = etime(t22,t11);

figure,plot(1:length(para_our.reg_Err),para_our.reg_Err,'r.-');hold on;
plot(1:length(para_our.reg_Diff),para_our.reg_Diff,'b.-');
set(gca,'YLim',[0,1]);

%% majiayi tip
tic;
[Transform1, C] = GLS_ma(X,Y,100,1,1);
toc;
reg_err_gls = measurement_XY(Transform1.Y,Y,Y,opt_gm.GT);

%% CPD
opt.method='rigid';
opt.corresp = 1;      % compute correspondence vector at the end of registration (not being estimated by default)
% XT = Y;
opt.viz = 1;
opt.inl_num = LX_sample;
[Transform,C_cpd]=cpd_register(Y,X,opt);

figure,plot(Transform.Y(:,1),Transform.Y(:,2),'bo',Y(:,1),Y(:,2),'r.');

reg_err_cpd = measurement_XY(Transform.Y,Y,Y,opt_gm.GT);


%% GAGM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_gagm = GAGM_similar_2d(X,Y,opt_reg,opt_gm);

%% IFFPS
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_ipfps = IPFPS_similar_2d(X,Y,opt_reg,opt_gm);

%% PSM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_psm = PSM_similar_2d(X,Y,opt_reg,opt_gm);

%% RRWM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_rrwm = RRWM_similar_2d(X,Y,opt_reg,opt_gm);

%% FGMD
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_fgmd = FGMD_similar_2d(X,Y,opt_reg,opt_gm);

%% MPM
opt_gm.GT = GT;
opt_gm.full_or_del = 'del';
opt_gm.topk = nF;

opt_reg.reg_normalize = 1;
opt_reg.display = 1;
opt_reg.write = 0;
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_mpm = MPM_similar_2d(X,Y,opt_reg,opt_gm);


%% FRGM
opt_gm = set_opt_FRGM(10,1,1);
opt_gm.GT = GT;
opt_gm.topk = nF;
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
opt_reg.rota  = 1;
opt_reg.reg_maxiter = 20;
opt_reg.unary_terms = [zeros(1,3),ones(1,50)];

para_frgm = FRGM_similar_2d_fast(X,Y,opt_reg,opt_gm);
