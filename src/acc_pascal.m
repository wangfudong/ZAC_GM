function [acc1,acc2] = acc_pascal(map,gt,nF)

[m,n] = size(map);
[r,c] = find(map==1);
[r1,c1]=find(gt==1);
right = 0;
right1 = 0;
for i=1:length(r1)
    right = right+ sum(ismember([r,c],[r1(i),c1(i)],'rows'));
    right1 = right1+ sum(ismember([r1,c1],[r(i),c(i)],'rows'));
end
acc1 = right/nF;
acc2 = right/min(m,n);