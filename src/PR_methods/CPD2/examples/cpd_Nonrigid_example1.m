%  Nonrigid Example1. Coherent Point Drift (CPD).
%  Registration of 2D fish point sets without noise and outliers.

%clear all; close all; clc;

% dataset = ['house';'hotel';'Carss';'Motor'];
% filenum = 2;
% if filenum <= 2
%     inds = [1,100];
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


load('cpd_data2D_fish.mat');
theta = -0.*pi;
Rota = [cos(theta),sin(theta);
    -sin(theta),cos(theta)];
Y = Y*Rota-1;
% order = randperm(length(X(:,1)));
% X = X(order,:);
opt.method='nonrigid';
opt.corresp = 1;
[Transform, C]=cpd_register(X,Y, opt);
axis equal
figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
