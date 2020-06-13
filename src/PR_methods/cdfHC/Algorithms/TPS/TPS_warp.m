%% Perform thin-plate spline warping
%% Input:
%%       landmarks:   source 2D pts stored in nx2 matrix.  
%%       parameters:  parameters in nx2 matrix where first three rows are
%%       affine parameters corresponding to <1,x,y>
%% Output:
%%       warped_pts:  target 2D pts in nx2 matrix
%%       energy:      bending energy

function [warped_pts, energy] = TPS_warp(landmarks, parameters)

[n,d] = size(landmarks);
[B,lambda] = compute_basis(landmarks);

warped_pts = B*parameters;
energy = trace(parameters(d+2:n,:)'*diag(lambda)*parameters(d+2:n,:));

