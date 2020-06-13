function Map = kLAP_Hun(M,k)
[m,n] = size(M);

if k < m
    M_tmp = zeros(m,m+n-k);
    M_tmp(1:m,1:n) = M;
    diagvec = zeros(1,m-k);
    for j = 1:m-k
        diagvec(j) = min(M(j,:))-1;
    end
    rightup = max(M(:))*m*n*ones(m-k,m-k);
    for i = 1:m-k
        rightup(i,i) = diagvec(i);
    end
    rightdown = repmat(diagvec,k,1);
    M_tmp(1:m-k,n+1:n+m-k) = rightup;
    M_tmp(m-k+1:m,n+1:n+m-k) = rightdown;
    
%     if exist('matchpairs')==2% if matlab version is 2019 or later
%         M0 = max(M_tmp(:));
%         Map1 = matchpairs(M_tmp,M0);
%         Map = zeros(size(M_tmp));
%         for i = 1:m
%             Map(Map1(i,1),Map1(i,2)) = 1;
%         end
%     else
%         Map = asgHun(-M_tmp);
%     end

    Map = asgHun(-M_tmp);

%     M_tmp = M_tmp/max(M_tmp(:)); % lapjv runs faster than asgHun with m,n >= 100
%     rowsol = lapjv(M_tmp,1);
%     Map = zeros(size(M_tmp));
%     for i = 1:m
%         Map(i,rowsol(i)) = 1;
%     end
    
    Map = Map(1:m,1:n);
else
    
    Map = asgHun(-M);
    
end