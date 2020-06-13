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

function [energys,after_tps] = obj_HC2_TPS_energy(param, init_affine, basis, kronBasis, kernel, scene, scale, alpha, beta,d, display)

display_it = display;
global colors;

npts=max(size(scene));

for i=1:npts
    [nLb{i},nb{i}] = size(basis{i});
end

stag=1;
for i=1:npts
    n=nb{i};
    transform_param{i}=param(stag:d*n+stag-1);
    stag=d*n+stag;
end

for i=1:npts
    n=nb{i};
    if isempty(init_affine{i})
        affine_param{i} = reshape(transform_param{i}(1:d*(d+1)),d,d+1);
        affine_param{i} = affine_param{i}';
        tps_param{i} = reshape(transform_param{i}(d*(d+1)+1:d*n),d,n-d-1);
        tps_param{i} = tps_param{i}';
    else
        tps_param{i} = reshape(transform_param{i}(1:d*n-d*(d+1)),d,n-d-1);
        tps_param{i} = tps_param{i}';
        affine_param{i} = reshape(init_affine{i},d,d+1);
        affine_param{i} = affine_param{i}';
    end
    after_tps{i} = basis{i}*[affine_param{i};tps_param{i}];
    bending(i)= trace(tps_param{i}'*kernel{i}*tps_param{i});
end

% after_tps{npts}=scene{npts};
% this is for biased registration

[ after_tps ] = TranslateToRPlus(after_tps,10*ones(1,d));

[energy] = HC2_energy(after_tps);

energys = alpha*sum(energy) + beta*sum(sum(bending));

%for display
if ( display_it == 1 )
    subplot(1,2,2);hold off;
    plot(after_tps{1}(:,1),after_tps{1}(:,2),'ro','linewidth',1.5,'markersize',6);
    hold on,plot(after_tps{2}(:,1),after_tps{2}(:,2),'b.','markersize',10);hold off;
    %     DisplayPoints_HC(after_tps,colors);
    %     title('After Registration','fontsize',16);
    drawnow;
end






