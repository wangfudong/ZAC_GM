% Compute the TPS parameters from two equally sized sets of 2D-points
% http://elonen.iki.fi/code/tpsdemo/index.html
% Based mostly on "Approximation Methods for Thin Plate Spline Mappings and Principal Warps" by Gianluca Donato and Serge Belongie, 2002.
function [A,w,energy] = compute_TPS(S, T)

[n,d] = size(S);

P = ones(n,d+1);
P(:,2:d+1) = S(:,1:d);

K = compute_K(S);

v = zeros(n+d+1,d);
v(1:n,:) =  T(:,1:d);

L = zeros(n+d+1,n+d+1);
L(1:n,1:n) = K;
L(1:n, n+1:n+d+1) = P;
L(n+1:n+d+1, 1:n) = P';

aw = L \ v;
w = aw(1:n,:);
A = aw(n+1:n+d+1,:);

energy = trace(w'*K*w);
