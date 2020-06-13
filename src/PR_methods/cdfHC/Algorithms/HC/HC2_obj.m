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

function [energy,grad]=HC2_obj(S,basis,alpha)

if nargin<3
    alpha=1/max(size(S));
end

npts=max(size(S));
d=min(size(S{1}));

energy1=0;
energy2=0;

for i=1:npts
    energy2=energy2+alpha*IntP2(S{i},S{i});
end

for i=1:npts
    for j=1:npts
        energy1=energy1+alpha^2*IntP2(S{i},S{j});
    end
end

energy=(-energy1+energy2);

for i=1:npts
    gradtemp{i}=zeros(max(size(S{i})),d);
end

for i=1:npts
    for j=1:npts
        if j~=i
           gradtemp{i}=gradtemp{i}+2*alpha^2*GradP2(S{i},S{j},2);      
        end
    end
    
    grad{i}=-gradtemp{i}+(alpha-alpha^2)*GradP2(S{i},S{i},1);
    
end








