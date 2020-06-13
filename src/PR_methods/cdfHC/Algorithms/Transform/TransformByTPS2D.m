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

% todo: add one more option, with-affine or not     
function [after_tps, bending] = TransformByTPS2D(before, param, with_affine)

[n,dim] = size(before);

if with_affine > 0
   param_length = 2*n+6;
   affine_param = param(1:6);
   tps_param = param(7:2*n+6);
else
   param_length = 2*n;
   affine_param = [1,0,0,1,0,0];
   tps_param = param(1:2*n);
end

if length(param) < param_length
    disp('At least 2n+6 parameters required for 2D TPS transform');
    return;
end

if (dim<2)
    disp('Input point sets should have dimensionality >=2');
    return;
end

if with_affine > 0
	A = reshape(affine_param,2,3);  
	%  [ a11 a12  | tx ]
	%  [ a21 a22  | ty ]
	%   note: order of reshape is column first
	%   param = [a11 a21 a12 a22 tx ty]
	after_affine = (A * [before(:,1:2)'; ones(1,n)])';
else
        after_affine = before;	
end
% if the input has more than 2 cols, the 3rd col will be treated as orientation
% information in degree 
% if dim>2
%      after(:,3) = before(:,3) + d_theta*180/pi;
% end        

W = reshape(tps_param, n, 2);
K = ones(n);
for i=1:n
    for j=1:n
        r = norm(before(i,1:2) - before(j,1:2));
        if (r==0)
            K(i,j) = 0;
        else
            K(i,j) =   r*r*log(r)/(8*pi);
        end
    end
end

after_tps = after_affine + K*W;
bending = abs(trace(W'*K*W));
