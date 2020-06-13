function opt = set_option(regular,inl_num,GT)

switch regular
    case 'w/o'% without regularization
        opt.unary = 1;
        opt.alpha1 = 1*1;
        opt.alpha2 = 1*1;
        opt.maxiter = 100;
        opt.miniter = 3;
        opt.active = 1;
        opt.q_norm = 2;
        opt.mass = inl_num;
        opt.GT = GT;
        opt.a = 0.5;
        opt.b = 0.5;
        opt.obj = 1;
        opt.remove_time = 3;
        opt.del_or_full = 'full';
        opt.output = 0;
        
    case 'w'
        opt.unary = 1;
        opt.alpha1 = 1*1;
        opt.alpha2 = 1*1;
        opt.alpha_r = 1;
        opt.maxiter = 100;
        opt.miniter = 3;
        opt.active = 1;
        opt.q_norm = 2;
        opt.q_norm_r = 2;
        opt.mass_r = inl_num;
        opt.GT = GT;
        opt.a = 0.5;
        opt.b = 0.5;
        opt.obj = 2;
        opt.remove_time = 3;
        opt.del_or_full = 'full';
        opt.output = 0;
        opt.r_fixed = 0;
              
end
