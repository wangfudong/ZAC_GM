function [Transform, C] = GLS_ma(X,Y,maxiter,display,verbose)

opt.outliers = 0.9;
opt.viz = display;
opt.t = 0.3;
opt.sparse = 0;
opt.nsc = 5;
opt.verbose = verbose;
opt.max_it = maxiter;
opt.normalize = 1;
opt.beta = 2;
opt.lambda = 3;
opt.tol = 1e-10;

[Transform, C]=prgls_register(Y, X, opt);
V = Transform.Y;

% if display >0
% figure,cpd_plot_iter(X, Y); axis off; title('Before');
% figure,cpd_plot_iter(V, Y); axis off; title('After registering Y to X');
% end