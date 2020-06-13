
% choose number of points
n = 50;

% generate random 2D point set of size n
landmarks = rand(n,2);

% randomly generate TPS parameters where first 3 rows are affine parameters
p = rand(n,2);

% warp
[warped_pts, energy1] = TPS_warp(landmarks, p);


% try to recover the parameters
% [A,w,energy2] = compute_TPS(landmarks, warped_pts);
% energy1
% energy2

% landmarks = rand(n,3);
% p = rand(n,3);
% [warped_pts, energy1] = TPS_warp(landmarks, p);
% [A,w,energy2] = compute_TPS(landmarks, warped_pts);
% energy1
% energy2
% A
% p'

m = 12;
n = 10;
landmarks = rand(m,3);
ctrl_pts = rand(n,3);

p = rand(n,3);
[warped_pts, energy1] = TPS_warp2(landmarks, ctrl_pts, p);
[A,w,energy2] = compute_TPS(landmarks, warped_pts);
energy1
energy2
A-p(1:4,:)