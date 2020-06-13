function [dis,u,v]=sinkhorn_NEG(a,b,K,U,stop,p_norm,tolerance,maxiter,verbose)
% compute entropy regularized OT distance and lower bounds on the EMD.
% Inputs:
%  a: m*1 or m*k, marginal measure  
%     - m*1: 1-vs-k model,dis=[d(a,b_1),...,d(a,b_k)]
%     - m*k: k times 1-vs-1 model ,dis=[d(a_1,b_1),...,d(a_k,b_k)]
%  b: n*k, marginal measure
%  M: m*n, distance between site_a and site_b
%  stop: stopping criteria, marginal-difference or convergence-decrese 
%  p_norm: used in stopping creteria
%  tolerence: stop creteria
%  iter: max iteration 
%  verbose: display the iteration 
% Outputs:
%  dis: entropy regularized OT distance
%  u: m*k, left scaling 
%  v: n*k, right scaling
%% processing inputs
if (nargin<5 || isempty(stop))
    stop = 'MarginalDifference';
end
if (nargin<6 || isempty(p_norm))
    p_norm = inf;
end
if(nargin<7 || isempty(tolerance))
    tolerance = 1.0e-6;
end
if (nargin<8 || isempty(maxiter))
    maxiter = 5000;
end
if (nargin<9 || isempty(verbose))
    verbose = 0;
end

%% preprocessing: check the dimensions of a and b
if (size(a,2) > 1)   
    error('a should be n*1 vector');
end
if max(a(:))~=1 || max(b(:))~=1
    error('max(a) or max(b) should equals 1');
end
%% preprocessing: remove zero values of a(i) in 1-vs-n case
index = (a>0);
existzeros = false;
if ~all(index) %
    existzeros = true;
    a=a(index);
    K=K(index,:);
    U=U(index,:);
end
%ainvK=bsxfun(@rdivide,K,a); % K./[a,a,...,a]

%% iterator, initialization
cnt=0;
u=ones(size(a,1),size(b,2))/length(index);

if strcmp(stop,'RelativeDecrease')
    dis_hold=ones(1,size(b,2)); %initialization of vector of distances.
end

%% iterator
while cnt < maxiter
    u=a./(K*(b./(K'*u)));
    cnt=cnt+1;
    
    if (mod(cnt,20)==1 || cnt==maxiter)
        v=b./(K'*u);
        u=a./(K*v);
        
        switch stop
            case 'RelativeDecrease'
                dis=sum(u.*(U*v));
                Criterion=norm(dis./dis_hold-1,p_norm);
                if Criterion<tolerance || isnan(Criterion),
                    break;
                end
                dis_hold=dis;
            case 'MarginalDifference',
                Criterion=norm(sum(abs(v.*(K'*u)-b)),p_norm);
                if Criterion<tolerance || isnan(Criterion), % norm of all || . ||_1 differences between the marginal of the current solution with the actual marginals.
                    break;
                end
            otherwise
                error('Stopping Criterion not recognized');
        end
        cnt=cnt+1;
        if verbose>0,
            disp(['Iteration :',num2str(cnt),' Criterion: ',num2str(Criterion)]);
        end
        if any(isnan(Criterion)), % stop all computation if a computation of one of the pairs goes wrong.
            error('NaN values have appeared during the fixed point iteration. This problem appears because of insufficient machine precision when processing computations with a regularization value of lambda that is too high. Try again with a reduced regularization parameter lambda or with a thresholded metric matrix M.');
        end
    end
end
  disp(['Iteration :',num2str(cnt),' Criterion: ',num2str(Criterion)]);  
 if strcmp(stop,'MarginalDifference'), % if we have been watching marginal differences, we need to compute the vector of distances.
    dis=sum(u.*(U*v));
 end

if nargout>2 && existzeros, % user wants scalings. We might have to arficially add zeros again in bins of a that were suppressed.
    uu=u;
    u=zeros(length(index),size(b,2));
    u(I,:)=uu;
end   
















