function dist = dist_desc(dx,dy)

[n, d]=size(dx); [m, d]=size(dy);

dist=repmat(dx,[1 1 m])-permute(repmat(dy,[1 1 n]),[3 2 1]);
dist=squeeze(sum(dist.^2,2));
dist = sqrt(dist)/d;