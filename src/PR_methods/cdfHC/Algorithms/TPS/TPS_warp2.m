%% Perform thin-plate spline warping
%% Input:
%%       landmarks:   source 2D pts stored in nx2 matrix.  
%%       parameters:  parameters in nx2 matrix where first three rows are
%%       affine parameters corresponding to <1,x,y>
%% Output:
%%       warped_pts:  target 2D pts in nx2 matrix
%%       energy:      bending energy

function [warped_pts, energy] = TPS_warp2(landmarks, ctrl_pts, parameters)


[m,d] = size(landmarks);
[n,d] = size(ctrl_pts);
[K,U] = compute_K(ctrl_pts,landmarks);

Pm = [ones(m,1) landmarks];
Pn = [ones(n,1) ctrl_pts];

PP = null(Pn'); B = [Pm U*PP]; 
warped_pts = B*parameters;
energy = trace(parameters(d+2:n,:)'*PP'*K*PP*parameters(d+2:n,:));
%% transformed = KKPP*rand(n,d);
%% bending = trace(w'*PP'*K*PP*w)

