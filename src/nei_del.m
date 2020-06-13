function nei = nei_del(X)
LX = length(X(:,1));
[ntriX] = delaunay(X(:,1),X(:,2));
nei = zeros(LX,LX);
for i = 1:length(ntriX(:,1))
    tri3 = sort(ntriX(i,:),'ascend');
    nei(tri3(1),tri3(2)) = 1;
    nei(tri3(1),tri3(3)) = 1;
    nei(tri3(2),tri3(3)) = 1;
end
