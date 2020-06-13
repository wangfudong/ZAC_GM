
% Generate model points

M=1500;

M1=round(M/3);
M2=M1;
M3=M-M1-M2;

model=zeros(3,M);
model(1:2,1:M1)=rand(2,M1);
model([1,3],(M1+1):(M1+M2))=rand(2,M2);
model(2,(M1+1):(M1+M2))=1;
model(2:3,(M1+M2+1):M)=rand(2,M3);
model(1,(M1+M2+1):M)=1;

% Generate data points

N=500;

N1=round(N/3);
N2=N1;
N3=N-N1-N2;

data=zeros(3,N);
data(1:2,1:N1)=rand(2,N1);
data([1,3],(N1+1):(N1+N2))=rand(2,N2);
data(2,(N1+1):(N1+N2))=1;
data(2:3,(N1+N2+1):N)=rand(2,N3);
data(1,(N1+N2+1):N)=1;

% Transform data points to their start positions

v1=0.6*(2*rand-1); v2=0.6*(2*rand-1); v3=0.6*(2*rand-1);
R1=[1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
R2=[cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
R3=[cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];

R=R3*R2*R1;

data=R*data;
data(1,:)=data(1,:)+0.2*randn;
data(2,:)=data(2,:)+0.2*randn;
data(3,:)=data(3,:)+0.2*randn;

% A plot. Model points and data points in start positions

figure(1)
plot3(model(1,:),model(2,:),model(3,:),'r.',data(1,:),data(2,:),data(3,:),'b.'), hold on, axis equal
plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
title('Original data points (blue) and model points (red)')

% Running the ICP-algorithm. Least squares criterion

[RotMat,TransVec,dataOut]=icp2(model,data);

% Reference:
%
% Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
% Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2


% A plot. Model points and data points in transformed positions

figure(2)
plot3(model(1,:),model(2,:),model(3,:),'r.',dataOut(1,:),dataOut(2,:),dataOut(3,:),'g.'), hold on, axis equal
plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
title('Transformed data points (green) and model points (red)')

