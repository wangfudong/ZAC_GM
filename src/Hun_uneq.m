function Map = Hun_uneq(M,mass)

[m,n] = size(M);
HP = asgHun(-M);

[row,col] = find(HP==1);
val = M(HP==1);

[~,ind] = sort(val,'ascend');
ind_r = row(ind(1:mass));
ind_c = col(ind(1:mass));
Map = zeros(m,n);
for j = 1:length(ind_r)
    Map(ind_r(j),ind_c(j)) = 1;
end


