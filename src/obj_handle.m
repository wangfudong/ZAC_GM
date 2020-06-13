function F_handle = obj_handle(M0,W1,W2,A,B,option)

unary = option.unary;
alpha1 = option.alpha1;
alpha2 = option.alpha2;
q2 = option.q_norm;

switch option.obj
        
    case 1
        F_handle = @(P)(unary*sum(sum(M0.*P)) +...
            alpha1*sum(sum(W1.*abs(A-P*B*P').^q2)) +...
            alpha2*sum(sum(W2.*abs(P'*A*P-B).^q2)));        
        
    case 2
        alphar = option.alpha_r;
        qr = option.q_norm_r;
        r = option.mass_r;
        
        F_handle = @(P)(unary*sum(sum(M0.*P)) +...
            alpha1*sum(sum(W1.*abs(A-P*B*P').^q2)) +...
            alpha2*sum(sum(W2.*abs(P'*A*P-B).^q2)))+...
            alphar*abs(sum(sum(P))-r)^qr;

end

