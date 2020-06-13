% MATLAB demo code of Max-Pooling Matching CVPR 2014
%
% M. Cho, J. Sun, O. Duchenne, J. Ponce
% Finding Matches in a Haystack: A Max-Pooling Strategy for Graph Matching in the Presence of Outliers 
% Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (2014) 
% http://www.di.ens.fr/willow/research/maxpoolingmatching/
%
% Please cite our work if you find this code useful in your research. 
%
% written by Minsu Cho, Inria - WILLOW / Ecole Normale Superieure 
% http://www.di.ens.fr/~mcho/

clear all; close all; clc;
disp('************************ Point Matching Test ************************');disp(' ');

%% Settings Evaluations
setPath; setPointMatching; setMethods;

%%
plotSet.lineWidth = 3; % Line width
plotSet.markerSize = 10; % Marker Size
plotSet.fontSize = 15; % Font Size
%plotSet.font = '\fontname{times new roman}'; % Font default
plotSet.font = '\fontname{Arial}'; % Font default

%% Test Methods
Accuracy = zeros(length(settings{Con}{4}), length(methods), Set.nTest);
MatchScore = zeros(length(settings{Con}{4}), length(methods), Set.nTest);
Time = zeros(length(settings{Con}{4}), length(methods), Set.nTest);
MatchScoreRaw = zeros(length(settings{Con}{4}), length(methods), Set.nTest);
MatchScoreMP = zeros(length(settings{Con}{4}), length(methods), Set.nTest);

t_start = clock;
fprintf(['Experiment starts: ' num2str(t_start(4)) ':' num2str(t_start(5)) ':' num2str(round(t_start(6))) '\n']);
for kk = 1:Set.nTest, fprintf('Test: %d of %d ', kk, Set.nTest);
    for i = 1:length(settings{Con}{4})
        eval(['Set.' settings{Con}{3} '=' num2str(settings{Con}{4}(i)) ';']);
        problem = makePointMatchingProblem(Set);
        eval(['Set.' settings{Con}{3} '= settings{' num2str(Con) '}{4};']);
        for j = 1:length(methods)
            [Accuracy(i,j,kk) MatchScore(i,j,kk) Time(i,j,kk) tmpX tmpXraw MatchScoreRaw(i,j,kk) MatchScoreMP(i,j,kk)] ...
                = wrapper_GM(methods(j), problem);
        end
        fprintf('.');
    end
    clf; handleCount = 0; yData = mean(Accuracy(:,:,1:kk),3); yLabelText = 'Accuracy'; plotResults;
    str = ['Average accuracy from test 1 to test ' num2str(kk)]; title(str, 'FontSize', 16); drawnow;
    t_now = clock; elap = etime(t_now, t_start); t_end = add_time(t_start, elap/(kk)*(Set.nTest));
    fprintf(['  expected time to end' num2str(t_end(4)) ':' num2str(t_end(5)) ':' num2str(round(t_end(6))) '\n']);
end
clear i j k temp X Xbin ind p val str
close all

%% Plot Results
meanAccuracy = mean(Accuracy,3);
meanMatchScore = mean(MatchScore,3);
meanTime = mean(Time,3);
meanMatchScoreRaw = mean(MatchScoreRaw,3);
meanMatchScoreMP = mean(MatchScoreMP,3);

%%
handleCount = 0;
yData = meanAccuracy; yLabelText = 'accuracy'; plotResults;
yData = meanMatchScore; yLabelText = 'objective score'; plotResults;
yData = meanTime; yLabelText = 'time'; plotResults;