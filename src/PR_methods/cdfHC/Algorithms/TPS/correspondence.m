function [S,T]=correspondence(X,Y)
n=size(X,1);

for i=1:n
    [x,y]=mindist(X(i,1),X(i,2),Y);
    T(i,1)=x;
    T(i,2)=y;
end
S=X;

    
    
function [x,y]=mindist(x0,y0,Y);
n=size(Y,1);
mind=dist(x0,y0,Y(1,1),Y(1,2));
for k=1:n
    d=dist(x0,y0,Y(k,1),Y(k,2));
    if mind>=d
        mind=d;
        count=k;
    end
end
x=Y(count,1);
y=Y(count,2);

    
function d=dist(x0,y0,x1,y1);
d=sqrt((x0-x1)^2+(y0-y1)^2);