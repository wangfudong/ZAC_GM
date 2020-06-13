function [distance, X] = GNCCP_A(Ct, Ag, Ah, para)
ng = size(Ag, 1);
nh = size(Ah, 1);
if ng <= nh
    if nargin > 3
        [distance, X] = GNCCP_GM(Ct, Ag, Ah, para);
    else
        [distance, X] = GNCCP_GM(Ct, Ag, Ah);
    end
else
    if nargin > 3
        [distance, X] = GNCCP_GM(Ct', Ah, Ag, para);
    else
        [distance, X] = GNCCP_GM(Ct', Ah, Ag);
    end
    X = X';
end
end