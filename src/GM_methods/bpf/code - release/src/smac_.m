function [D, P, K] = smac(Ag, Ah, K)
M = size(Ag, 1);    N = size(Ah, 1);
if M < N
   Ag_ins = zeros(N, N);
   Ag_ins(1:M, 1:M) = Ag;
   Ag = Ag_ins;
end
if nargin < 3
    K = get_K(Ag, Ah);
end
% refer to balanced graph matching supplementary materials
constraint_col = kron(eye(N, N), ones(1, N));
constraint_row = kron(ones(1, N), eye(N, N));
C = [constraint_col; constraint_row];
b = ones(2*N, 1);
Ceq = [eye(2*N-1) zeros(2*N-1, 1)] * (C - b * C(2*N, :));
PC = eye(N*N) - Ceq'*inv(Ceq*Ceq' + 0.0001 * eye(size(Ceq, 1)))*Ceq;
Kac = PC * K * PC;

[p, eig] = eigs(Kac, 1);
P = reshape(abs(p), N, N);
P = P';
P = KM(P);
P = P(1:M,:);

D = sum(sum((Ag - P * Ah * P').*(Ag - P * Ah * P')));

function K = get_K(Ag, Ah)
siz_g = size(Ag, 1);
siz_h = size(Ah, 1);
K = sparse(siz_g * siz_h, siz_g * siz_h);
for i = 1 : siz_g
    for j = 1 : siz_g
        K(((i - 1) * siz_h + 1) : i * siz_h , ((j - 1) * siz_h + 1) : j * siz_h) = -(Ag(i, j) - Ah).^2;
    end
end
minK = min(min(K));
K = K - minK;

function P = permutation(vec, m, n)
P = zeros(m, n);
for i = 1 : m
   P(i, vec(i)) = 1; 
end   

