% Example 1. Rigid CPD point-set registration. No options are set, so the
% default ones are used. 2D fish point-set.
clear all; close all; clc;

load cpd_data2D_fish; 
% delete some points and add outliers.
X = cpdfish{1, 1};[X,~] = normalize_point(X,1);
%Y = cpdfish{1, 2};[Y,~] = normalize_point(Y,1);
% Add a random rotation and scaling
R=cpd_R(0.2*pi);
s=1.5;
Y=s*X*R';

X = [X(10:end,:);0.3*rand(20,2)];
Y = [Y(1:(end-10),:);0.1*rand(20,2)];

opt.fgt = 1;
Transform=cpd_register(X,Y,opt);

% Initial point-sets
figure,cpd_plot_iter(X, Y); title('Before');

% Registered point-sets
figure,cpd_plot_iter(X, Transform.Y);  title('After');

% Rotation and scaling errors after the registration
E_R=norm(R-Transform.R)
E_s=norm(s-Transform.s)
