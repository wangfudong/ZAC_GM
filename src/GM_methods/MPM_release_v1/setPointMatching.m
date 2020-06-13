%% Configuration Settings
nSet=0; %% Display (bool) / Display Name / Variable Name / Variable Value, Vector, String
nSet=nSet+1; settings{nSet} = {0, '# of tests', 'nTest', 10};
switch 1
    case 1 % Outlier Test
        nSet=nSet+1; settings{nSet} = {1, '# of inliers \itn_{\rmin}\rm', 'nInlier', 10};
        nSet=nSet+1; settings{nSet} = {1, '# of outliers \itn_{\rmout}\rm', 'nOutlier', 0:10:100};
        nSet=nSet+1; settings{nSet} = {1, 'deformation noise \it\sigma\rm', 'deformation', 0.03};
    case 2 % Deformation Test
        nSet=nSet+1; settings{nSet} = {1, '# of inliers \itn_{\rmin}\rm', 'nInlier', 20};
        nSet=nSet+1; settings{nSet} = {1, '# of outliers \itn_{\rmout}\rm', 'nOutlier', 0};
        nSet=nSet+1; settings{nSet} = {1, 'deformation noise \it\sigma\rm', 'deformation', 0:0.02:0.2};
    
    otherwise
end

nSet=nSet+1; settings{nSet} = {0, 'Point Distribution', 'typeDistribution', 'normal'}; % normal / uniform
nSet=nSet+1; settings{nSet} = {0, 'Transformation Scale', 'transScale', 1}; % scale change
nSet=nSet+1; settings{nSet} = {0, 'Transformation Rotation', 'transRotate', 0}; % rotation change
nSet=nSet+1; settings{nSet} = {0, 'Scale for 2nd order', 'scale_2D', 0.5}; % Used in 2D
nSet=nSet+1; settings{nSet} = {0, 'Scale for 3rd order', 'scale_3D', 0.5}; % Used in 3D
nSet=nSet+1; settings{nSet} = {0, 'Permute Nodes', 'bPermute', 1}; % boolean
nSet=nSet+1; settings{nSet} = {0, 'Outliers Both Side', 'bOutBoth', 0}; % boolean
nSet=nSet+1; settings{nSet} = {0, 'Use Displacement', 'bDisplacement', 0}; % boolean

%% Evaluate Settings
nFix = 0; nCon = 0;
for n = 1:nSet
    if isscalar(settings{n}{4})
        eval(['Set.' settings{n}{3} '=' num2str(settings{n}{4}) ';']);    
        if settings{n}{1}, nFix=nFix+1; Fix(nFix) = n; end
    else
        eval(['Set.' settings{n}{3} '= settings{n}{4};']);    
        if settings{n}{1}, nCon=nCon+1; Con(nCon) = n; end
    end
end
clear n nFix nCon nSet

disp('* Check experiment settings *'); disp(Set);