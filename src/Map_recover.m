function Map = Map_recover(map,m,n,xind,yind)
Map = zeros(m,n);
[m1,n1] = size(map);
for i=1:m1
    for j = 1:n1
        Map(xind(i),yind(j)) = map(i,j);
    end
end