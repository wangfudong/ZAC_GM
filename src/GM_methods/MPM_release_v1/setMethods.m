%% Methods & Settings
% Script for setting algorithms to run

% You can add an algorithm following the script below
%nMethods = 1;
%methods(nMethods).fhandle = @fhandle;                         % Function of the algorithm
%methods(nMethods).variable = {'var1', 'var2', 'var3'};        % Input variables that the algorithm requires
%methods(nMethods).param = {'name1', 'val1', 'name2', 'val2'}; % Default parameter values
%methods(nMethods).strName = 'algorithm name';                 % Algorithm name tag
%methods(nMethods).color = 'color';                            % Color for plots
%methods(nMethods).lineStyle = 'line style';                   % Line style for plots
%methods(nMethods).marker = 'marker';                          % Marker for plots

nMethods = 0;
if 1
    nMethods = nMethods + 1;
    methods(nMethods).fhandle = @MPM;
    methods(nMethods).variable = {'affinityMatrix', 'group1', 'group2'};
    methods(nMethods).param = {};
    methods(nMethods).strName = 'MPM';
    methods(nMethods).color = 'r';
    methods(nMethods).lineStyle = '-';
    methods(nMethods).marker = 'o';
end
if 1
    nMethods = nMethods + 1;
    methods(nMethods).fhandle = @RRWM;
    methods(nMethods).variable = {'affinityMatrix', 'group1', 'group2'};
    methods(nMethods).param = {};
    methods(nMethods).strName = 'RRWM';
    methods(nMethods).color = 'b';
    methods(nMethods).lineStyle = '-';
    methods(nMethods).marker = 's';
end
if 1
    nMethods = nMethods + 1;
    methods(nMethods).fhandle = @ipfp_gm;
    methods(nMethods).variable = {'affinityMatrix', 'L12'};
    methods(nMethods).strName = 'IPFP';
    methods(nMethods).param = {};
    methods(nMethods).color = 'm';
    methods(nMethods).lineStyle = '-';
    methods(nMethods).marker = 'x';
end
if 1
    nMethods = nMethods + 1;
    methods(nMethods).fhandle = @SM;
    methods(nMethods).variable = {'affinityMatrix'};
    methods(nMethods).param = {};
    methods(nMethods).strName = 'SM';
    methods(nMethods).color = 'k';
    methods(nMethods).lineStyle = '--';
    methods(nMethods).marker = 'x';
end

%% Show the algorithms to run
disp('* Algorithms to run *');
for k = 1:nMethods, disp([methods(k).strName ' : @' func2str(methods(k).fhandle)]); end; disp(' ')
clear k