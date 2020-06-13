function [E, G] = costfun(param, X, Y, K, U, lambda, sigma2, is_grad)

[N, D] = size(X);
M = size(U, 2);
C = reshape(param, [M D]);

E = 1/( 2^D * pi^(D/2) * sigma2^(D/4) );
E = E + lambda/2 * trace(C'*K*C);

V = U*C - Y;
a = -2 / N / (2*pi*sigma2)^(D/2);
F = exp(-sum(V.^2, 2) / (2*sigma2));
E = E + a * sum(F);

%%
G = [];
if is_grad
    G = - a * U' * ( V .* repmat(F, [1 D]) / sigma2 ) + lambda * K * C;
    G = G(:);
end

% [N, D] = size(X);
% C = reshape(param, [N D]);
% 
% E = 1/( 2^D * pi^(D/2) * sigma2^(D/4) );
% E = E + lambda/2 * trace(C'*K*C);
% 
% a = -2 / N / (2*pi*sigma2)^(D/2);
% 
% tmp = zeros(1,N);
% for n=1:N
%     tmp(n) = exp(-sum((K(n,:)*C-Y(n,:)).^2)/(2*sigma2));
% end
% E = E + a * sum(tmp);
% 
% %%
% G = [];
% if is_grad
%     G = zeros(N,2);
%     for n=1:N
%        G = G - tmp(n)*K(n,:)'*(K(n,:)*C-Y(n,:))/sigma2; 
%     end
%     G = a * G + lambda * K * C;
%     G = G(:);
% end

