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

function gm=gradymin(ti,tj,beta)
% if nargin<3
%     beta=50;
% end
% gm=exp(-beta*ti)/(exp(-beta*ti)+exp(-beta*tj));

if ti<=tj
    gm=1;
else gm=0;
end