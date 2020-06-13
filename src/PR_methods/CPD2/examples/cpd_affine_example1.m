% Example 1. Affine CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear all; close all; clc;

% dataset = ['house';'hotel';'Carss';'Motor'];
% filenum = 2;
% if filenum <= 2
%     inds = [1,85];
% else
%     inds = 15;
% end
% [Pts,Ims,nF] = fileload(filenum,inds);
% X = Pts{1,1};
% Y = Pts{1,2};
% I1 = Ims{1,1};
% I2 = Ims{1,2};
% if filenum <= 2
%     LX = 30;
%     LY = 30;
%     order = randperm(nF);
% order = order(1:LX);
% %order = [27,26,10,17,2,3,24,12,23,25,6,29,30,5,28,7,9,20,22,16];
% remove = 0;
% else
%     outliers = 20;
%     LX = nF;
%     LY = min(LX + outliers,length(Y(:,1)));
%     order = 1:LX;
%     remove = 1;
% end
% X = X(order,:);

load cpd_data2D_fish; Y=X;

%%Add a random affine transformation
  B=eye(2)+0.5*abs(randn(2,2));
X=X*B'-3;


opt.method='affine';
opt.corresp=1;      % compute correspondence vector at the end of registration (not being estimated by default)

Transform=cpd_register(X,Y,opt);

% Initial point-sets
figure,cpd_plot_iter(X, Y); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X, Transform.Y);  title('After');
