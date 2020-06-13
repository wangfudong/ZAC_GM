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

function [ s ] = TranslateToRPlus(s,region)

n = max(size(s));
a = [s{1}];
for i = 2:n
    a = [a;s{i}];
end
for i = 1:n
    s{i} = s{i}-repmat(min(a)-region,size(s{i},1),1);
end