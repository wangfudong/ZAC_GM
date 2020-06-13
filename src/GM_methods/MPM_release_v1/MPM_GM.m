function asgMPM = MPM_GM(X,Y,order,K_del,scale_2D)


nP1 = length(order);
nP2 = length(Y(:,1));
P1 = X(order,:);
P2 = Y;
E12 = ones(nP1,nP2);
n12 = nnz(E12);

[L12(:,1), L12(:,2)] = find(E12);
[group1, group2] = make_group(L12);

E1 = ones(nP1); E2 = ones(nP2);
[L1(:,1), L1(:,2)] = find(E1);
[L2(:,1) ,L2(:,2)] = find(E2);
G1 = P1(L1(:,1),:)-P1(L1(:,2),:);
G2 = P2(L2(:,1),:)-P2(L2(:,2),:);


G1 = sqrt(G1(:,1).^2+G1(:,2).^2);
G2 = sqrt(G2(:,1).^2+G2(:,2).^2);
G1 = reshape(G1, [nP1 nP1]);
G2 = reshape(G2, [nP2 nP2]);
M = (repmat(G1, nP2, nP2)-kron(G2,ones(nP1))).^2;


M = exp(-M./scale_2D);
M(1:(n12+1):end)=0;

gth = zeros(nP1,nP2);
for i = 1:nP1
    gth(i,order(i))=1;
end
GT.bool = gth(:);

problem.nP1 = nP1;
problem.nP2 = nP2;
problem.P1 = P1;
problem.P2 = P2;
problem.L12 = L12;
problem.E12 = E12;

if nargin < 4 || isempty(K_del)
    problem.affinityMatrix = M;
else
    problem.affinityMatrix = full(K_del);
end
%problem.affinityMatrix = M;

problem.group1 = group1;
problem.group2 = group2;

problem.GTbool = GT.bool;

% setMethods;

methodtest.fhandle= @MPM;
methodtest.variable= {'affinityMatrix'  'group1'  'group2'};
methodtest.param= {};
methodtest.strName= 'MPM';
methodtest.color='r';
methodtest.lineStyle= '-';
methodtest.marker='o';

[~,score,time,X,X_con] = wrapper_GM(methodtest, problem);

asgMPM.X_con = X_con;
asgMPM.obj = score;
asgMPM.tim = time;
asgMPM.X = reshape(X,[nP1,nP2]);
[Maprate,~] = ratemap(asgMPM.X,order,1:nP2,0);
asgMPM.acc = Maprate;
