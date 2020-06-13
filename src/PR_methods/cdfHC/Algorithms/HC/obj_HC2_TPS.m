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

function [energys,gradient] = obj_HC2_TPS(param, init_affine, basis, kronBasis, kernel, scene, scale, alpha, beta,d,display)
% obj_HC2_TPS(x,     affine,  TPS_basis, kronBasis, TPS_kernel, S, scale, alpha, beta, d,display_it)
global colors;

display_it = display;
npts=length(scene);

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

[ after_tps ] = TranslateToRPlus(after_tps,[1,1]);

[energy,grad]=HC2_obj(after_tps,basis);

energys = alpha*sum(energy) + beta*sum(sum(bending));

%gradient computation
for i=1:npts
    nL=nLb{i};
    n=nb{i};
    grad{i} = alpha*basis{i}'*grad{i}; %kron(eye(d),basis);
    grad{i}(d+2:n,:) = grad{i}(d+2:n,:) + 2*beta*kernel{i}*tps_param{i};
    if isempty(init_affine{i})
        grad{i} = grad{i}';
        grad{i} = reshape(grad{i},1,d*n);
    else
        grad{i}(1:d+1,:) = [];
        grad{i} = grad{i}';
        grad{i} = reshape(grad{i},1,d*(n-d-1));
    end
end
gradient=[grad{1}];
for i=2:npts
    gradient=[gradient,grad{i}];
end
gradient=gradient';

%for display
if ( display_it == 1 )
    subplot(1,2,2);
    plot(after_tps{1}(:,1),after_tps{1}(:,2),'ro','linewidth',1.5,'markersize',6);
    hold on,plot(after_tps{2}(:,1),after_tps{2}(:,2),'b.','markersize',10);hold off;
    %DisplayPoints_HC(after_tps,colors);
    %title('After Registration','fontsize',16);
    drawnow;
end


