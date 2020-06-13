%% Group-wise Point-set registration using a novel CDF-based Havrda-Charvat Divergence (CDF-HC) Demo (version 20111209):
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

%% Please cite this paper if you use any pieces of the code.
%
% Reference: Ting Chen, Baba C. Vemuri, Anand Rangarajan and Stephan J. Eisenschenk,
%            Group-wise Point-set registration using a novel CDF-based Havrda-Charvat Divergence.
%            In IJCV : International Journal of Computer Vision, 86(1):111-124, January, 2010.
%% To test this script, run the any of the following
% load cdfhc_data2D_beijing
% HC2Reg_TPS(beijing,1,200);
% 
load cpd_data2D_fish
cpdfish{1,2} = cpdfish{1,2}-repmat([1,1],size(cpdfish{1,2},1),1);
[T,En] = HC2Reg_TPS(cpdfish,1,200);
% 
% load cdfhc_data2D_CC
% HC2Reg_TPS(CC7,1,200);

