
function [B,lambda] = compute_basis(landmarks)

[n,d] = size(landmarks);
K = compute_K(landmarks);


P = [ones(n,1) landmarks];
Q = inv(K)*P;
%p_basis = orth(P);
q_basis = orth(Q);
%[U,S,V] = svd(Q);
%Q_perp_basis = orth(U(:,4:n));
A = eye(n) - q_basis*q_basis';
G = A*inv(K)*A;
[V,D] = eig(G);
g_basis = V(:,1:n-d-1);
B = [P,g_basis];
lambda = diag(D(1:n-d-1,1:n-d-1));

%% null(P') perp Im(P);
%% PP = null(P'); KPP = K*PP; KKPP = [P K*PP]; 
%% transformed = KKPP*rand(n,d);
%% bending = trace(w'*PP'*K*PP*w)