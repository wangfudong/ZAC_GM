function X = GAM_outlier_1p4(M,group1,group2,XGT)

a = 1.4;
nIn = sum(XGT(:));
n1 = size(group1,2);
n2 = size(group2,2);
assert(n1==n2);
nOut = n1 - nIn;

E12_slack = ones(nIn+nOut*2);
[L12_slack(:,1), L12_slack(:,2)] = find(E12_slack);
[group1_slack, group2_slack] = make_group12(L12_slack);
Mslack = a*ones((nIn+nOut*2)^2);
id = ismember(L12_slack(:,1),1:n1) & ismember(L12_slack(:,2),1:n2);
Mslack(id,id) = M;

% Xraw_slack = RRWM_upload(double(Mslack), group1_slack, group2_slack);
Xraw_slack = GAM(double(Mslack),L12_slack,group1_slack,group2_slack);
Xslack = zeros(size(E12_slack)); Xslack(find(E12_slack)) = Xraw_slack;
Xslack = discretisationMatching_hungarian(Xslack,E12_slack); Xslack = Xslack(find(E12_slack));

Xslack = reshape(Xslack,n1+nOut,n2+nOut);
X = Xslack(1:n1,1:n2);
X = X(:);

% only for bPermute==0
% fprintf('sum(X) = %f, sum(Xtrue) = %f\n',sum(X),sum(sum(Xslack(1:nIn,1:nIn))));

end