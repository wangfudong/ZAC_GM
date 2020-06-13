function [Y]=tpswarp(S,T)

x=S(:,1);
y=S(:,2);
n=size(x,1);
one=ones(n,1);
P=[one,x,y];
Pt=P';
K=zeros(n,n);
for i=1:n
    for j=1:n
        K(i,j)= U(sqrt((x(i)-x(j))^2+(y(i)-y(j))^2));
    end
end


A=[K,P;Pt,zeros(3,3)];

V=[T;zeros(3,2)];
p=inv(A)*V;

X=A*p;
Y=X(1:n,:);