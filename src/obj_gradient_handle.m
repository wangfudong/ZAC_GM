function gF = obj_gradient_handle(M,W1,W2,A,B,option)

unary = option.unary;
alpha1 = option.alpha1;
alpha2 = option.alpha2;
[m,n] = size(M);
q2 = option.q_norm;


switch option.obj
    case 1 
        P = [];
        q1 = option.q_prob;
        small_eps = 1.0e-5;
        L1 = sum(P,2).^q1*(sum(P,2)').^q1;
        L2 = (sum(P,1)').^q1*sum(P,1).^q1;
        W1_tmp1 = W1.*abs(A-P*B*P').^q2;
        W1_tmp2 = -q2*W1.*L1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P');
        W2_tmp1 = W2.*abs(P'*A*P-B).^q2;
        W2_tmp2 = q2*W2.*L2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B);
        gF1 = zeros(size(P));
        gF2 = gF1;
        if q1 >= 1
            for i = 1:m
                gF1(i,:) = q1*(sum(P(i,:))^(q1-1))*repmat(2*W1_tmp1(i,:)*sum(P,2),1,n) + 2*W1_tmp2(i,:)*(P*B');
            end
        elseif q1 < 1
            for i = 1:m
                gF1(i,:) = q1*(max(sum(P(i,:)),small_eps)^(q1-1))*repmat(2*W1_tmp1(i,:)*sum(P,2),1,n) + 2*W1_tmp2(i,:)*(P*B');
            end
        end
        if q1 >= 1
            for i = 1:n
                gF2(:,i) = q1*(sum(P(:,i)).^(q1-1))*repmat(2*W2_tmp1(i,:)*sum(P,1)',m,1) + 2*A*P*W2_tmp2(i,:)';
            end
        elseif q1 < 1
            for i = 1:n
                gF2(:,i) = q1*(max(sum(P(:,i)),small_eps).^(q1-1))*repmat(2*W2_tmp1(i,:)*sum(P,1)',m,1) + 2*A*P*W2_tmp2(i,:)';
            end
        end
        
    case 2
        P = [];
        a = option.a;
        b = option.b;
        L1 = (1+exp((a-sum(P,2))/b)).^(-1)*(1+exp((a-sum(P,2)')/b)).^(-1);
        L2 = (1+exp((a-sum(P,1)')/b)).^(-1)*(1+exp((a-sum(P,1))/b)).^(-1);
        
        W1_tmp1 = W1.*abs(A-P*B*P').^q2;
        W1_tmp2 = -q2*W1.*L1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P');
        W2_tmp1 = W2.*abs(P'*A*P-B).^q2;
        W2_tmp2 = q2*W2.*L2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B);
        gF1 = zeros(size(P));
        gF2 = gF1;
        for i = 1:m
            gF1(i,:) = (1/b)*(1+exp((a-sum(P(i,:)))/b)).^(-2)*exp(a-sum(P(i,:)))*repmat(2*W1_tmp1(i,:)*(1+exp((a-sum(P,2))/b)).^(-1),1,n) + 2*W1_tmp2(i,:)*(P*B');
        end
        for i = 1:n
            gF2(:,i) = (1/b)*(1+exp((a-sum(P(:,i)))/b)).^(-2)*exp((a-sum(P(:,i)))/b)*repmat(2*W2_tmp1(i,:)*(1+exp((a-sum(P,1)')/b)).^(-1),m,1) + 2*A*P*W2_tmp2(i,:)';
        end
        
    case 3
        P = [];
        small_eps = 1.0e-20;
        L11 = max(sum(P,2),small_eps).^(-1);
        L22 = max(sum(P,1),small_eps).^(-1);
        L1 = exp(1-L11)*exp(1-L11)';
        L2 = exp(1-L22)'*exp(1-L22);
        
        W1_tmp1 = W1.*abs(A-P*B*P').^q2;
        W1_tmp2 = -q2*W1.*L1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P');
        W2_tmp1 = W2.*abs(P'*A*P-B).^q2;
        W2_tmp2 = q2*W2.*L2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B);
        gF1 = zeros(size(P));
        gF2 = gF1;
        for i = 1:m
            gF1(i,:) = (1-L11(i)*exp(1-L11(i))*L11(i).^2)*repmat(2*W1_tmp1(i,:)*exp(1-L11),1,n) + 2*W1_tmp2(i,:)*(P*B');
        end
        for i = 1:n
            gF2(:,i) = (1-L22(i)*exp(1-L22(i))*L22(i).^2)*repmat(2*W2_tmp1(i,:)*exp(1-L22)',m,1) + 2*A*P*W2_tmp2(i,:)';
        end
        
    case 4
        W1_tmp2 = @(P) (-q2*W1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P'));
        W2_tmp2 = @(P) (q2*W2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B));
        
        gF1 = @(P) (2*W1_tmp2(P)*(P*B'));
        gF2 = @(P) (2*A*P*W2_tmp2(P)');
        
%         gF1 = zeros(size(P));
%         gF2 = gF1;
%         for i = 1:m
%             gF1(i,:) = 2*W1_tmp2(i,:)*(P*B');
%         end
%         for i = 1:n
%             gF2(:,i) = 2*A*P*W2_tmp2(i,:)';
%         end
        
    case 5
        
        UV1_tmp = @(P) (sum(P,2)*(sum(P,2)'));
        UV2_tmp = @(P) ((sum(P,1)')*sum(P,1));
        W1_tmp2 = @(P) (-q2*W1.*UV1_tmp(P).*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P'));
        W2_tmp2 = @(P) (q2*W2.*UV2_tmp(P).*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B));
        
        gF1 = @(P) (2*W1_tmp2(P)*(P*B'));
        gF2 = @(P) (2*A*P*W2_tmp2(P)');
        
        %         gF1 = zeros(size(P));
        %         gF2 = gF1;
        %         for i = 1:m
        %             gF1(i,:) = 2*W1_tmp2(i,:)*(P*B');
        %         end
        %         for i = 1:n
        %             gF2(:,i) = 2*A*P*W2_tmp2(i,:)';
        %         end
        
    case 6
        W1_tmp2 = @(P) (-q2*W1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P'));
        W2_tmp2 = @(P) (q2*W2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B));
        
        gF1 = @(P) (2*W1_tmp2(P)*(P*B'));
        gF2 = @(P) (2*A*P*W2_tmp2(P)');
        
    case 7
        UV1_tmp = @(P) (sum(P,2)*(sum(P,2)'));
        UV2_tmp = @(P) ((sum(P,1)')*sum(P,1));
        W1_tmp2 = @(P) (-q2*W1.*UV1_tmp(P).*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P'));
        W2_tmp2 = @(P) (q2*W2.*UV2_tmp(P).*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B));
        
        gF1 = @(P) (2*W1_tmp2(P)*(P*B'));
        gF2 = @(P) (2*A*P*W2_tmp2(P)');
        
    case 8
        gF1 = @(P) (-2*(A*P*B' + A'*P*B));
        gF2 = gF1;
        
    case 9
        gF1 = @(P) (-2*(A*P*B' + A'*P*B));
        gF2 = gF1;
end

gF = @(P) (unary*M + alpha1*gF1(P) + alpha2*gF2(P));

switch option.obj
    case 6
        gF = @(P) (gF(P) + option.alpha_r*option.q_norm_r*abs(sum(sum(P))-option.mass_r)^(option.q_norm_r-1)*sign(sum(sum(P))-option.mass_r));
    case 7
        gF = @(P) (gF(P) + option.alpha_r*option.q_norm_r*abs(sum(sum(P))-option.mass_r)^(option.q_norm_r-1)*sign(sum(sum(P))-option.mass_r));
    case 9
        gF = @(P) (gF(P) + option.alpha_r*option.q_norm_r*abs(sum(sum(P))-option.mass_r)^(option.q_norm_r-1)*sign(sum(sum(P))-option.mass_r));
end




