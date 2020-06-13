function  [X, Y, normal] =norm2(x,y)
% NORM2 nomalizes the data so that it has zero means and unit covariance

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

x = double(x);
y = double(y);

n=size(x,1);
m=size(y,1);

normal.xm=mean(x);
normal.ym=mean(y);

x=x-repmat(normal.xm,n,1);
y=y-repmat(normal.ym,m,1);

normal.xscale=sqrt(sum(sum(x.^2,2))/n);
normal.yscale=sqrt(sum(sum(y.^2,2))/m);

X=x/normal.xscale;
Y=y/normal.yscale;