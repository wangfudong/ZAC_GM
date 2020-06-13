%% HC point set registration with TPS 
function [T,energys]=HC2Reg_TPS(S, display, MaxIter, scale, alpha, beta, init_p, init_affine)

% S       : point sets to be registered
% display : 0--not display, 1--display, 1 is the default value
% MaxIter : number of iterations
% scale   : scale of the point sets, set to be 1 as default 
% alpha, beta : regularization parameters, alpha=1, beta=0.00001 as default
% init_p, init_affine : parameter initialization 
 
%% main
global colors;

d = min(size(S{1}));
npts = length(S); 
colors = rand(npts,3);
%colors = [1,0,0;0,0,1];

if nargin<5
    scale = 1; 
    alpha = 1;
    beta = 0.00001;
end

if nargin<3
    MaxIter = 80;
end

if nargin < 2
    display = 1;
end

% move point sets to R+ region
[ S ] = TranslateToRPlus(S,[1,1]);


for i=1:npts
    n{i} = max(size(S{i})); % n{i}: number of points in S{i} 
end

if nargin < 8
    for i = 1:npts
        init_affine{i} = [0, 0, 1, 0, 0, 1]';   %[tx,ty,a11,a12,a21,a22];   
    end
end

for i=1:npts
      x0{i}= [init_affine{i}; zeros(d*n{i}-d*(d+1),1)];  % 0 initialization for all parameters 
      affine{i} = [ ];
end 

% for display
display_it = display;
if ( display_it == 1 )
%     subplot(1,2,1);hold off;
%     plot(S{1}(:,1),S{1}(:,2),'r.')
%     DisplayPoints_HC(S,colors); 
%     title('Initial Configuration','fontsize',16);
%     drawnow;
    
    hh = figure;subplot(1,2,1);
   set(hh,'position',[500 450 1000 300]);
    plot(S{1}(:,1),S{1}(:,2),'r.',S{2}(:,1),S{2}(:,2),'b.','markersize',20);
    %hold on,plot(Y(LX+1:end,1),Y(LX+1:end,2),'g.','markersize',10);
    drawnow;
end  

% compute TPS basis, kernel
for i=1:npts
    [m{i},d] = size(S{i});
    [n{i},d] = size(S{i});
    [K{i},U{i}] = compute_K(S{i},S{i});
    Pm{i} = [ones(m{i},1) S{i}];
    Pn{i} = [ones(n{i},1) S{i}];
    PP{i} = null(Pn{i}'); 
    TPS_basis{i} = [Pm{i} U{i}*PP{i}]; 
    TPS_kernel{i} = PP{i}'*K{i}*PP{i};
    kronBasis{i} = kron(eye(d),TPS_basis{i});
    disp('TPS:basis,kernel computed'); 
end

%options = optimset('GradObj','on','Display', 'none','LargeScale','off', 'MaxIter', MaxIter,  'TolFun',1e-012, 'TolX',1e-015);
% %%Display the iterations
options = optimset('GradObj','on','Display', 'none','LargeScale','off', 'MaxIter', MaxIter,  'TolFun',1e-012, 'TolX',1e-015);

% re-arrange initial param
X0=[x0{1}];
for i=2:npts 
    X0=[X0;x0{i}];
end

if nargin >= 7
    X0=init_p;
end

t = cputime;
[param, fval, exit_flag, output] = fminunc(@(x)obj_HC2_TPS(x, affine, TPS_basis, kronBasis, TPS_kernel, S, scale, alpha, beta, d,display_it), X0,  options);

timecost = cputime - t

[energys,T] = obj_HC2_TPS_energy(param, affine, TPS_basis, kronBasis, TPS_kernel, S, scale, alpha, beta, d,display_it);

