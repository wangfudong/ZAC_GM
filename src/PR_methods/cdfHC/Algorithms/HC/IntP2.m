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

function v=IntP2(Pts1,Pts2)

n1=max(size(Pts1));
n2=max(size(Pts2));
d=min(size(Pts1));
v=0;

for i=1:n1
    for j=1:n2
       %v=v+analymin(Pts1(i,1),Pts2(j,1))*analymin(Pts1(i,2),Pts2(j,2)); 
       v=v+min(Pts1(i,1),Pts2(j,1))*min(Pts1(i,2),Pts2(j,2));
    end
end

v=1/(n1*n2)*v;