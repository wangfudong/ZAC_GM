function F = obj_value(M,P,W1,W2,A,B,option)

unary = option.unary;
alpha1 = option.alpha1;
alpha2 = option.alpha2;
q2 = option.q_norm;

switch option.obj
    case 1
        q1 = option.q_prob;
        
        L1 = sum(P,2).^q1*(sum(P,2)').^q1;
        L2 = (sum(P,1)').^q1*sum(P,1).^q1;
        F = unary*sum(sum(M.*P)) + alpha1*sum(sum(W1.*L1.*abs(A-P*B*P').^q2)) + alpha2*sum(sum(W2.*L2.*abs(P'*A*P-B).^q2));
        
    case 2
        a = option.a;
        b = option.b;
        L1 = ((1+exp((a-sum(P,2))/b)).^(-1)*(1+exp((a-sum(P,2)')/b)).^(-1));
        L2 = ((1+exp((a-sum(P,1)')/b)).^(-1)*(1+exp((a-sum(P,1))/b)).^(-1));
        F = (unary*sum(sum(M.*P)) +...
            alpha1*sum(sum(W1.*L1.*abs(A-P*B*P').^q2)) +...
            alpha2*sum(sum(W2.*L2.*abs(P'*A*P-B).^q2)));
        
    case 3
        L1 = (exp(1-sum(P,2).^(-1))*exp(1-sum(P,2).^(-1))');
        L2 = (exp(1-sum(P,1).^(-1))'*exp(1-sum(P,1).^(-1)));
        F = (unary*sum(sum(M.*P)) +...
            alpha1*sum(sum(W1.*L1.*abs(A-P*B*P').^q2)) +...
            alpha2*sum(sum(W2.*L2.*abs(P'*A*P-B).^q2)));
    case 4
        F = (unary*sum(sum(M.*P)) +...
            alpha1*sum(sum(W1.*abs(A-P*B*P').^q2)) +...
            alpha2*sum(sum(W2.*abs(P'*A*P-B).^q2)));
end




