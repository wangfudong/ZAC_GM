function [ viewInfo ] = extract_localfeatures_wrapper( filePathName, tmpPathName, varargin )
set_param_GM;
% load input data
tmpPath = [ tmpPathName '_view.mat' ];
if exist(tmpPath)
    load( tmpPath, 'viewInfo');
else
    fprintf('- extracting local features from %s ...\n', filePathName);
    viewInfo = extract_localfeatures_mcho( filePathName, fparam );
    if ~isempty(tmpPathName)
        save( tmpPath, 'viewInfo');
    end
end