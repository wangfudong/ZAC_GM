% -------------------------------------------------------------------
% Copyright (C) 2006-2011 Ting Chen, Baba C. Vemuri and Anand Rangarajan
% 
% Authors: Ting Chen
% Date:    12/09/2011
% 
% Contact Information:
%
% Ting Chen, tichen@cise.ufl.edu
% 
% Terms:	  
% 
% The source code (M-files) are provided under the
% terms of the GNU General Public License with an explicit
% clause permitting the execution of the M-files from within
% a MATLAB environment. See the LICENSE file for details.
% ------------------------------------------------------------------

function DisplayPoints_HC(S,colors)

n = length(S);

plot(S{1}(:,1),S{1}(:,2),'o',...
               'MarkerEdgeColor',colors(1,:),...
               'MarkerFaceColor', colors(1,:),...
               'MarkerSize',6);
for i =2:n
    hold on; plot(S{i}(:,1),S{i}(:,2),'o',...
               'MarkerEdgeColor',colors(i,:),...
               'MarkerFaceColor', colors(i,:),...
               'MarkerSize',6);
end

axis equal;
axis tight;
set(gca,'FontSize',12);
axis off;



