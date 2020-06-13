function A = DelaunayTriGraph(kps)
% Build Delaunay Trianglelation graph.
%
% Input
%   kps   -  n * 2 matrix of coordinates of keypoints
%  
% Output
%   A     - weighted adjcent matrix of graph
%
% History
%   create   -  Tao WANG (twang@bjtu.edu.cn), 12-28-2014

    %% max distance between keypoints
    maxd = 0.0;
    for i = 1:size(kps, 1)
        for k = 1:size(kps, 1)
            dist = getKeypointDist(kps, i, k);
            if (dist > maxd)
                maxd = dist;
            end
        end
    end
%     maxd = 1.5 * maxd;
    
    %% build Delaunay Triangle graph
    x = kps(:, 1);
    y = kps(:, 2);
    n = size(kps, 1);
    TRI = delaunay(x, y);
    
    %% build adjcent matrix 
    A = zeros(n, n);
    for i = 1:size(TRI, 1)
        A(TRI(i,1), TRI(i, 2)) = 1.0 - getKeypointDist(kps, TRI(i,1), TRI(i, 2)) / maxd;
        A(TRI(i,1), TRI(i, 3)) = 1.0 - getKeypointDist(kps, TRI(i,1), TRI(i, 3)) / maxd;
        A(TRI(i,2), TRI(i, 3)) = 1.0 - getKeypointDist(kps, TRI(i,2), TRI(i, 3)) / maxd;
%         A(TRI(i,1), TRI(i, 2)) = getKeypointDist(kps, TRI(i,1), TRI(i, 2)) / maxd;
%         A(TRI(i,1), TRI(i, 3)) = getKeypointDist(kps, TRI(i,1), TRI(i, 3)) / maxd;
%         A(TRI(i,2), TRI(i, 3)) = getKeypointDist(kps, TRI(i,2), TRI(i, 3)) / maxd;
% 
        A(TRI(i,2), TRI(i, 1)) = A(TRI(i,1), TRI(i, 2));
        A(TRI(i,3), TRI(i, 1)) = A(TRI(i,1), TRI(i, 3));
        A(TRI(i,3), TRI(i, 2)) = A(TRI(i,2), TRI(i, 3));
     end
end

function dist = getKeypointDist(kps, i, k)
    x = kps(i,1) - kps(k,1);
    y = kps(i,2) - kps(k,2);
    dist = sqrt(x * x + y * y);
end