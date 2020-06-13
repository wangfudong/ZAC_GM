function cdata = conWillowGphs(img1, anno1, img2, anno2, outliers)

cdata = ExtractFeature(img1, anno1, img2, anno2, outliers);

cdata = BuildAffinityMatrix(cdata);

%cdata = BuildFgmUMatrices(cdata);

cdata = BuildFgmDMatrices(cdata);


%assert(abs(sum(cdata.KU(:)-cdata.KD(:))) < eps);

%KK = cdata.K - cdata.KD;

end

function cdata = ExtractFeature(img1, anno1, img2, anno2, nOutliers)
set_param_GM;

viewInfo1 = extract_localfeatures_mcho( img1, fparam, 'verbose', true );
load(anno1, 'pts_coord');%% load the annotated data    
gts1 = size(pts_coord,2);

% annotations 1
nPts1 = size(viewInfo1.frame, 2);
idx_all1 = 1:nPts1;
idx_sel1 = [];
for i=1:gts1  %size(pts_coord,2)
    pt1 = pts_coord(1:2,i);
    %idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'Radius', distClose);
    idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'NumberOfNeighbours', 1);
    idx_sel1 = [ idx_sel1, idx1 ];
    
    idx_all1(i) = idx1;
    idx_all1(idx1) = i;
end
% add outliers
idx_sel1 = [idx_sel1, idx_all1(i + randperm(nPts1 - i, nOutliers))];
cdata.Pts1 = viewInfo1.frame(1:2, idx_sel1)';
n1 = size(cdata.Pts1, 1);

viewInfo1.type = viewInfo1.type( idx_sel1 );
viewInfo1.frame = viewInfo1.frame( :, idx_sel1);
viewInfo1.desc = viewInfo1.desc( :, idx_sel1);


viewInfo2 = extract_localfeatures_mcho( img2, fparam, 'verbose', true );
load(anno2, 'pts_coord');%% load the annotated data    
gts2 = size(pts_coord, 2);

% annotations 2
nPts2 = size(viewInfo2.frame, 2);
idx_all2 = 1:nPts2;
idx_sel2 = [];
for i=1:gts2   %size(pts_coord,2)
    pt2 = pts_coord(1:2,i);
    %idx2 = nearestneighbour(pt2, viewInfo2.frame(1:2,:), 'Radius', distClose);
    idx2 = nearestneighbour(pt2, viewInfo2.frame(1:2,:), 'NumberOfNeighbours', 1);
    idx_sel2 = [ idx_sel2 idx2 ];
    
    idx_all2(i) = idx2;
    idx_all2(idx2) = i;    
end
%add outliers
idx_sel2 = [idx_sel2, idx_all2(i + randperm(nPts2 - i, nOutliers))];
cdata.Pts2 = viewInfo2.frame(1:2, idx_sel2)';
n2 = size(cdata.Pts2, 1);

viewInfo2.type = viewInfo2.type( idx_sel2 );
viewInfo2.frame = viewInfo2.frame( :, idx_sel2);
viewInfo2.desc = viewInfo2.desc( :, idx_sel2);

cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
cdata.view(1) = viewInfo1;  cdata.view(2) = viewInfo2;

assert(gts1 == gts2);

cdata.nGts = gts1;
cdata.X_GT = zeros(n1, n2);
cdata.X_GT(1:gts1, 1:gts2) = eye(gts1);
end

function cdata = BuildAffinityMatrix(cdata)
%% compute affinity between all possible matches
[ matchInfo ]= make_initialmatches_mcho2( cdata.view, cdata.mparam, 'verbose', true );
cand_matchlist = matchInfo.match;
cand_matchdist = matchInfo.dist;


[ uniq_feat1, tmp, new_feat1 ] = unique(cand_matchlist(1,:));    
[ uniq_feat2, tmp, new_feat2 ] = unique(cand_matchlist(2,:));
matchlist_tmp = [ new_feat1'; new_feat2' ];

% edgeAttr1 = computeEdgeAttr( cdata.view(1).frame(:,uniq_feat1) );
% edgeAttr2 = computeEdgeAttr( cdata.view(2).frame(:,uniq_feat2) );
% featBin1 = makeFeatBin(edgeAttr1);
% featBin2 = makeFeatBin(edgeAttr2);
% [ eSimVal ] = computeDotProdSimilarity_sym( int32(matchlist_tmp), featBin1, featBin2); % symmetric affinities


[ cdata.group1, cdata.group2 ] = make_group12(matchlist_tmp(1:2,:));
%[ cdata.L12, cdata.E12, cdata.nP1, cdata.nP2 ] = updateL12E12( matchlist_tmp' );

cdata = rmfield(cdata, 'view');

cdata.K = zeros(size(cand_matchdist, 2));
cdata.KD = zeros(size(cand_matchdist, 2));

%unaryWeight = 10.0*ones(size(cand_matchdist));
unaryWeight = 1.0*ones(size(cand_matchdist));
vSimVal = max( 0, 0.8 - cand_matchdist );
nodeAff = unaryWeight(matchlist_tmp(1,:)) .* vSimVal;
try
    cdata.K(1:(size(cdata.K,1)+1):end)= nodeAff;
    cdata.KD(1:(size(cdata.KD,1)+1):end)= nodeAff;
catch
    error('\n fail to assign node affinity!');
end

cdata.matchInfo = matchInfo.match;

end

function cdata = BuildFgmUMatrices(cdata)
cdata.gphs{1}.Pts = cdata.Pts1;
cdata.gphs{2}.Pts = cdata.Pts2;

%% build graphs
n1 = size(cdata.gphs{1}.Pts, 1);
n2 = size(cdata.gphs{2}.Pts, 1);

A1 = DelaunayTriGraph(cdata.gphs{1}.Pts);
A2 = DelaunayTriGraph(cdata.gphs{2}.Pts);

%% edges
[p1, q1, ~] = find(tril(A1));
[p2, q2, ~] = find(tril(A2));
cdata.gphs{1}.Eg = [p1',q1'; q1',p1'];
cdata.gphs{2}.Eg = [p2',q2'; q2',p2'];

% matrices G and H for undirect FGM-U
m1 = size(p1,1);    
m2 = size(p2,1);
cdata.gphs{1}.G = zeros(n1, m1);
cdata.gphs{2}.G = zeros(n2, m2);
for i = 1:m1
    cdata.gphs{1}.G(p1(i), i) = 1;    
    cdata.gphs{1}.G(q1(i), i) = 1;
end
for i = 1:m2
    cdata.gphs{2}.G(p2(i), i) = 1;
    cdata.gphs{2}.G(q2(i), i) = 1;
end
cdata.gphs{1}.H = [cdata.gphs{1}.G eye(n1)];
cdata.gphs{2}.H = [cdata.gphs{2}.G eye(n2)];


%% node affinity
cdata.KP = zeros(n1, n2);
nodeAff = cdata.K(1:(size(cdata.K,1)+1):end);
for i = 1 : n1*n2
    cdata.KP(cdata.matchInfo(1,i), cdata.matchInfo(2,i)) = nodeAff(i);
end

%% edge affinity
match_ind = zeros(n1, n2);
for i = 1:size(cdata.matchInfo, 2)
    match_ind(cdata.matchInfo(1,i), cdata.matchInfo(2,i)) = i;
end


S1 = nanstd(cdata.gphs{1}.Pts(:, 1:2));
S2 = nanstd(cdata.gphs{2}.Pts(:, 1:2));
cdata.KQ = zeros(m1, m2);

for i = 1:m1
    u1 = cdata.gphs{1}.Eg(1,i);
    v1 = cdata.gphs{1}.Eg(2,i);

    cor_1u = cdata.gphs{1}.Pts(u1, 1:2);
    cor_1v = cdata.gphs{1}.Pts(v1, 1:2);
    y_dist = abs(cor_1u(1) - cor_1v(1)) / S1(1);
    x_dist = abs(cor_1u(2) - cor_1v(2)) / S1(2);
    dst1 = sqrt(y_dist * y_dist + x_dist * x_dist);
    ang1 = atan(y_dist / (x_dist + eps));

    for k = 1:m2
        u2 = cdata.gphs{2}.Eg(1,k);
        v2 = cdata.gphs{2}.Eg(2,k);
        
        cor_2u = cdata.gphs{2}.Pts(u2, 1:2);
        cor_2v = cdata.gphs{2}.Pts(v2, 1:2);
        y_dist = abs(cor_2u(1) - cor_2v(1)) / S2(1);
        x_dist = abs(cor_2u(2) - cor_2v(2)) / S2(2);
        dst2 = sqrt(y_dist * y_dist + x_dist * x_dist);
        ang2 = atan(y_dist / (x_dist + eps));
        
        Dst = abs(dst1 - dst2);
        Ang = abs(abs(ang1) - abs(ang2));
        
        cdata.KQ(i,k) = exp(-(Dst + Ang) / 2);
        
        ind_u1u2 = match_ind(u1, u2);      % u1 --> u2
        ind_u1v2 = match_ind(u1, v2);      % u1 --> v2 
        ind_v1v2 = match_ind(v1, v2);      % v1 --> v2
        ind_v1u2 = match_ind(v1, u2);      % v1 --> u2
        
        % symmetric for undirected graphs
        cdata.K(ind_u1u2, ind_v1v2) = cdata.KQ(i,k);
        cdata.K(ind_v1v2, ind_u1u2) = cdata.KQ(i,k);
     %   cdata.K(ind_u1v2, ind_v1u2) = cdata.KQ(i,k);
     %   cdata.K(ind_v1u2, ind_u1v2) = cdata.KQ(i,k);
        
    end
end

cdata.K = sparse(cdata.K);

%cdata.x_gt = (matchInfo.match(1,:) == matchInfo.match(2,:)) & (matchInfo.match(1,:) <= cdata.nGts);
%cdata.x_gt = cdata.x_gt';

end

function cdata = BuildFgmDMatrices(cdata)
cdata.gphDs{1}.Pts = cdata.Pts1;
cdata.gphDs{2}.Pts = cdata.Pts2;

%% build graphs
n1 = size(cdata.gphDs{1}.Pts, 1);
n2 = size(cdata.gphDs{2}.Pts, 1);

A1 = DelaunayTriGraph(cdata.gphDs{1}.Pts);
A2 = DelaunayTriGraph(cdata.gphDs{2}.Pts);

%% edges
[p1, q1, ~] = find(tril(A1));
[p2, q2, ~] = find(tril(A2));
cdata.gphDs{1}.Eg = [p1',q1'; q1',p1'];
cdata.gphDs{2}.Eg = [p2',q2'; q2',p2'];

% matrices G and H for direct FGM-D 
p1 = cdata.gphDs{1}.Eg(1,:)';
q1 = cdata.gphDs{1}.Eg(2,:)';
p2 = cdata.gphDs{2}.Eg(1,:)';
q2 = cdata.gphDs{2}.Eg(2,:)';
m1 = size(p1,1);    
m2 = size(p2,1);
cdata.gphDs{1}.G = zeros(n1, m1);
cdata.gphDs{1}.H = zeros(n1, m1);
cdata.gphDs{2}.G = zeros(n2, m2);
cdata.gphDs{2}.H = zeros(n2, m2);
for i = 1:m1
    cdata.gphDs{1}.G(p1(i), i) = 1;    
    cdata.gphDs{1}.H(q1(i), i) = 1;
end
for i = 1:m2
    cdata.gphDs{2}.G(p2(i), i) = 1;
    cdata.gphDs{2}.H(q2(i), i) = 1;
end


%% node affinity
cdata.KP = zeros(n1, n2);
nodeAff = cdata.KD(1:(size(cdata.KD,1)+1):end);
for i = 1 : n1*n2
    cdata.KP(cdata.matchInfo(1,i), cdata.matchInfo(2,i)) = nodeAff(i);
end

%% edge affinity
match_ind = zeros(n1, n2);
for i = 1:size(cdata.matchInfo, 2)
    match_ind(cdata.matchInfo(1,i), cdata.matchInfo(2,i)) = i;
end


S1 = nanstd(cdata.gphDs{1}.Pts(:, 1:2));
S2 = nanstd(cdata.gphDs{2}.Pts(:, 1:2));
cdata.KQD = zeros(m1, m2);

for i = 1:m1
    u1 = cdata.gphDs{1}.Eg(1,i);
    v1 = cdata.gphDs{1}.Eg(2,i);

    cor_1u = cdata.gphDs{1}.Pts(u1, 1:2);
    cor_1v = cdata.gphDs{1}.Pts(v1, 1:2);
    y_dist = abs(cor_1u(1) - cor_1v(1)) / S1(1);
    x_dist = abs(cor_1u(2) - cor_1v(2)) / S1(2);
    dst1 = sqrt(y_dist * y_dist + x_dist * x_dist);
    ang1 = atan(y_dist / (x_dist + eps));

    for k = 1:m2
        u2 = cdata.gphDs{2}.Eg(1,k);
        v2 = cdata.gphDs{2}.Eg(2,k);
        
        if (i > m1 /2 && k <= m2 / 2) || (i <= m1 / 2 && k > m2 /2)
            continue;
        end
        
        cor_2u = cdata.gphDs{2}.Pts(u2, 1:2);
        cor_2v = cdata.gphDs{2}.Pts(v2, 1:2);
        y_dist = abs(cor_2u(1) - cor_2v(1)) / S2(1);
        x_dist = abs(cor_2u(2) - cor_2v(2)) / S2(2);
        dst2 = sqrt(y_dist * y_dist + x_dist * x_dist);
        ang2 = atan(y_dist / (x_dist + eps));
        
        Dst = abs(dst1 - dst2);
        Ang = abs(abs(ang1) - abs(ang2));
        
        cdata.KQD(i,k) = exp(-(Dst + Ang) / 2);
        
        ind_u1u2 = match_ind(u1, u2);      % u1 --> u2
    %    ind_u1v2 = match_ind(u1, v2);      % u1 --> v2 
        ind_v1v2 = match_ind(v1, v2);      % v1 --> v2
    %    ind_v1u2 = match_ind(v1, u2);      % v1 --> u2
   
        cdata.KD(ind_u1u2, ind_v1v2) = cdata.KQD(i,k);
    %    cdata.K(ind_v1v2, ind_u1u2) = cdata.KQD(i,k);
    %    cdata.K(ind_u1v2, ind_v1u2) = cdata.KQD(i,k);
    %    cdata.K(ind_v1u2, ind_u1v2) = cdata.KQD(i,k);
    
    end
end

cdata.KD = sparse(cdata.KD);

%cdata.x_gt = (matchInfo.match(1,:) == matchInfo.match(2,:)) & (matchInfo.match(1,:) <= cdata.nGts);
%cdata.x_gt = cdata.x_gt';

end

