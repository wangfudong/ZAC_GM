
% Generate model points

xvals=linspace(0,2*pi,100);
yvals=sin(xvals);

model=[xvals;yvals];

% Generate data points

xvals=linspace(0,2*pi,50);
yvals=sin(xvals);

data=[xvals;yvals];

% data = X';
% model = XT';

% Transform data points to their start positions

v1=0.6*(2*rand-1);
Rma=[cos(v1) -sin(v1);sin(v1) cos(v1)];

data=Rma*data;
data(1,:)=data(1,:)+2*randn;
data(2,:)=data(2,:)+2*randn;

% A plot. Model points and data points in start positions

figure(1)
plot(model(1,:),model(2,:),'r.',data(1,:),data(2,:),'b.'), axis equal

% Running the ICP-algorithm. Least squares criterion

[RotMat,TransVec,dataOut]=icp2(model,data);

% A plot. Model points and data points in transformed positions

figure(2)
plot(model(1,:),model(2,:),'r.',dataOut(1,:),dataOut(2,:),'b.'), axis equal

