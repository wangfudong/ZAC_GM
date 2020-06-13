function [R,s,t] = rigid_parameter(X0,Y0,P)
uX = sum(P,2)'*X0/sum(P(:));
uY = sum(P,1)*Y0/sum(P(:));
[m,n] = size(P);

X = X0-repmat(uX,m,1);
Y = Y0-repmat(uY,n,1);

M = X'*P*Y;
[U,~,V] = svd(M);
C = eye(2);
C(2,2) = det(U*V');
R = U*C*V';

s = trace(R'*X'*P*Y)/trace(X'*diag(sum(P,2))*X);

t = uY - s*uX*R;