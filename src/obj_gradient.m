function gF = obj_gradient(M,P,W1,W2,A,B,option)

unary = option.unary;
alpha1 = option.alpha1;
alpha2 = option.alpha2;
q2 = option.q_norm;


switch option.obj
    
    case 1
        W1_tmp2 = -q2*W1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P');
        W2_tmp2 = q2*W2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B);
        
        gF1 = 2*W1_tmp2*(P*B');
        gF2 = 2*A*P*W2_tmp2';
        
        
    case 2
        W1_tmp2 = -q2*W1.*abs(A-P*B*P').^(q2-1).*sign(A-P*B*P');
        W2_tmp2 = q2*W2.*abs(P'*A*P-B).^(q2-1).*sign(P'*A*P-B);
        
        gF1 = 2*W1_tmp2*(P*B');
        gF2 = 2*A*P*W2_tmp2';
        
        
end

gF = unary*M + alpha1*gF1 + alpha2*gF2;

if option.obj == 2
    gF = gF + option.alpha_r*option.q_norm_r*abs(sum(sum(P))-option.mass_r)^(option.q_norm_r-1)*sign(sum(sum(P))-option.mass_r);
    
end




