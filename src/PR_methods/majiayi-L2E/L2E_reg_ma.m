function [idt,V] = L2E_reg_ma(X,Y,maxiter,visualize,verbose)

normalize = 1;
%visualize = 1;

% configuration of SC
mean_dist_global=[]; % use [] to estimate scale from the data
nbins_theta=12;
nbins_r=5;
nsamp1=size(X,1);
nsamp2=size(Y,1);
ndum1=max(0, nsamp2-nsamp1);
ndum2=max(0, nsamp1-nsamp2);
eps_dum=0.15;
r_inner=1/8;
r_outer=2;
n_iter=maxiter;

if visualize
   hh = figure;subplot(1,2,1);
   set(hh,'position',[500 450 1000 300]);
    plot(X(:,1),X(:,2),'r.',Y(:,1),Y(:,2),'b.','markersize',20);
    %hold on,plot(Y(LX+1:end,1),Y(LX+1:end,2),'g.','markersize',10);
    drawnow;
end

% initialize transformed version of model pointset
Xk=X; 
% initialize counter
k=1;
s=1;
% out_vec_{1,2} are indicator vectors for keeping track of estimated
% outliers on each iteration
out_vec_1=zeros(1,nsamp1); 
out_vec_2=zeros(1,nsamp2);

while s
   disp(['iter=' int2str(k)])
   
   % compute shape contexts for (transformed) model
   %[BH1,mean_dist_1]=sc_compute(Xk',zeros(1,nsamp1),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_1);
   [BH1,mean_dist_1]=sc_compute_rota(Xk',zeros(1,nsamp1),mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec_1);

   % compute shape contexts for target, using the scale estimate from
   % the warped model
   % Note: this is necessary only because out_vec_2 can change on each
   % iteration, which affects the shape contexts.  Otherwise, Y does
   % not change.
   %[BH2,mean_dist_2]=sc_compute(Y',zeros(1,nsamp2),mean_dist_1,nbins_theta,nbins_r,r_inner,r_outer,out_vec_2);
   [BH2,mean_dist_2]=sc_compute_rota(Y',zeros(1,nsamp2),mean_dist_1,nbins_theta,nbins_r,r_inner,r_outer,out_vec_2);

   % compute pairwise cost between all shape contexts
   costmat=hist_cost_2(BH1,BH2);
   % pad the cost matrix with costs for dummies
   nptsd=nsamp1+ndum1;
   costmat2=eps_dum*ones(nptsd,nptsd);
   costmat2(1:nsamp1,1:nsamp2)=costmat;
   disp('running hungarian alg.')
   cvec=hungarian_ma(costmat2);
%   cvec=hungarian_fast(costmat2);
   [a,cvec2]=sort(cvec);
   disp('done.')
   
   X2 = Xk;
   Y2 = Y(cvec2(1:nsamp1),:);
   
   normal.xm=0; normal.ym=0;
   normal.xscale=1; normal.yscale=1;
   if normalize, [nX, nY, normal]=norm2(X2,Y2); end
   
   [idt, V] = L2E(nX, nY-nX, 0.5,verbose);
   if normalize, V=(V+nX)*normal.yscale+repmat(normal.ym,size(Y2,1),1); end 
   out_idx = setdiff(1:nsamp1, idt);
%    out_vec_1(out_idx) = 1;
%    out_vec_2(cvec2(out_idx)) = 1;
%    out_vec_2(cvec2(nsamp1+1:end)) = 1;
   
   if visualize
%       figure(2)
%       plot(V(:,1),V(:,2),'b+',Y(:,1),Y(:,2),'ro')
%       hold off
%       title([int2str(length(idt)) ' correspondences (warped X)'])
%       drawnow	    
        subplot(1,2,2); 
        plot(V(:,1),V(:,2),'ro','linewidth',1.5,'markersize',6);
        hold on,plot(Y(:,1),Y(:,2),'b.','markersize',10);
        hold off;
        drawnow;
   end
   
%    if visualize
%       figure(3)
%       plot(V(:,1),V(:,2),'b+',Y(:,1),Y(:,2),'ro')
%       hold on
%       h=plot([V(idt,1) Y2(idt,1)]',[V(idt,2) Y2(idt,2)]','k-');
%       hold off
%       title([int2str(length(idt)) ' correspondences (warped X)'])
%       drawnow	
%    end
%    
   % update Xk for the next iteration
   Xk = V;
   
   % stop early if shape context score is sufficiently low
   if k==n_iter
      s=0;
   else
      k=k+1;
   end
end
