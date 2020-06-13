function [sol, it, residual] = SSOR( mtx_A, rhs, init_vec, omega )

n       = size( mtx_A, 1 );
tol     = 10 * eps;
maxit   = 1000;
it      = 1;
flag    = 1;
tmp     = init_vec;
sol     = zeros(n,1);
%
% Take the strictly lower triangular matrix of matrix A
%
mtx_L   = tril(mtx_A, -1);
%
% Take the diagonal matrix of matrix A
%
diag_D  = diag(diag(mtx_A));
%
DomegaL = diag_D + omega * mtx_L;
DL      = ( 1 - omega ) * diag_D - omega * mtx_L;

%
% Compute M(omega)^(-1) * b
%
M_inv_b = DomegaL' \ ( diag_D * ( DomegaL \ rhs ) );
M_inv_b = ( omega * ( 2 - omega ) ) * M_inv_b;

while ( it <= maxit && flag == 1 )
    %
    % Compute T(omega) * tmp
    %
    sol = DomegaL' \ ( DL * ( DomegaL \ ( DL' * tmp ) ) );
    %
    % Compute new vector sol + M_inv_b
    sol = sol + M_inv_b;
    
    residual = norm(mtx_A * sol - rhs);
%     fprintf('residual(%4.0f) = %17.10e \n', it, residual);
    if ( residual <= tol || norm(sol-tmp)/norm(sol) <= tol )
        flag = 0;
    else
        tmp  = sol;
        it   = it + 1;
    end
end
