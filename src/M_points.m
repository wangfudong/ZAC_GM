function MP = M_points(X,Y,q)
if nargin < 3 || isempty(q)
    q = 0.5;
end
MP =  (bsxfun(@minus,X(:,1),Y(:,1)').^2 + bsxfun(@minus,X(:,2),Y(:,2)').^2).^q;
%MP = MP/max(MP(:));