%%sinkhorn algorithm

% tset among test_equal.m test_sinkhorn.m test_partial.m with equal constraint
da = 0.01:0.01:1;
db = 0.01:0.01:1;
mu_a = 0.3;delta_a=0.01;
mu_b = [0.1,0.3,0.7,0.9];delta_b=[0.01,0.01,0.05,0.08];

p = exp(-(da'-mu_a).^2/delta_a);p=p/sum(p);% m*1 列向量
b = zeros(length(db),length(mu_b));% n*k 
for i = 1:length(b(1,:))
    b(:,i) = exp(-(db'-mu_b(i)).^2/delta_b(i));
    b(:,i) = b(:,i)/sum(sum(b(:,i)));
end
q = b(:,3);
figure,plot(1:length(da),p,'r.');
figure,plot(1:length(db),q,'b.');


%% ground-distance
dis = abs(bsxfun(@minus,da',db)).^(2);%a在竖轴，b在横轴
dis = dis.^(1/1);
dis = dis/max(dis(:));
%dis = dis + 0.5;
%dis = dis./median(dis(:));
%figure,plot3(repmat(da',1,length(db)),repmat(db,length(da),1),dis,'r.');
figure,imagesc(dis);
%% shinkorn fixed point algorithm
lambda = 1/100;
K = exp(-1/lambda*dis);
figure,imagesc(K)
stop = 'RelativeDecrease';

tol = 1.0e-5;
maxiter = 5000;
verbose = 1;

tic;
[Map,Dis,u,v,iter,err_hold] = sinkhorn_OT(p,q,dis,lambda,tol,maxiter,verbose); % running with VERBOSE
toc;

%%
figure,imagesc(Map),title('map T');
figure,surf(Map),title('map T');

% Tr_norm=Map./repmat(p,1,length(q));
% Tc_norm=Map./repmat(q',length(p),1);
% figure,plot3(repmat(da',1,length(db)),repmat(db,length(da),1),Map,'r.');
% figure,imagesc(Tr_norm);title('Tr norm');
% figure,surf(Tr_norm);title('Tr norm');
% figure,imagesc(Tc_norm);title('Tc norm');
% figure,surf(Tc_norm),title('Tc norm');

tt = max(err_hold(:,1));
figure,
for i = 1:tt
    if err_hold(i,2) > 0
        subplot(1,2,1), plot(i,err_hold(i,2),'r.');hold on,
        subplot(1,2,2), plot(i,err_hold(i,3),'g.');hold on,
    end
end
