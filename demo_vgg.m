%% example of VGG dataset
datalist = {'graf';'bark';'bikes';'boat';'leuven';'trees';'ubc';'wall'};
graph_path = 'data\sift_graphs\';

data_id = 1;
dataname = datalist{data_id};

id1 = 1;
id2 = 3;

graph_mat = [graph_path dataname '_pair' num2str(id1) num2str(id2) '.mat'];
load(graph_mat);

XPoint = GRAPH.Xrest;
YPoint = GRAPH.Yrest;
descx = GRAPH.descrestX;
descy = GRAPH.descrestY;
Xmatch_result = GRAPH.Xmatch_result;
Ymatch_result = GRAPH.Ymatch_result;
maxsize1 = max(max(GRAPH.im1_size));
maxsize2 = max(max(GRAPH.im2_size));
H = GRAPH.X2Ymatrix;
pix_tol = GRAPH.pixel_th;
x_inlnum = Xmatch_result(1);
y_inlnum = Ymatch_result(1);

maxsize = max(maxsize1, maxsize2);
XX = XPoint/maxsize;
YY = YPoint/maxsize;

LX = length(XX(:,1));
LY = length(YY(:,2));
%% ZAC and ZACR
inl_rate = 0.5;
nF = round(inl_rate*min(LX,LY));

GT = [];
opt = set_option('w/o',nF,GT);
% opt = set_option('w',nF,GT);% for ZACR

opt.descx = descx;
opt.descy = descy;
opt.unary = 10;

tic;
[Map,output] = ZAC_vgg(XX,YY,opt);
% [Map,output] = ZAC_vgg_r(XX,YY,opt);%for ZACR
toc;

Map_cor = matrix2vec(asgHun(Map));
X_plot = XPoint(output.X_inl,:);
Y_plot = YPoint(output.Y_inl,:);

% [X2Y_correct,Y2X_correct,xre,xpre,yre,ypre] = match_result_sift(X_plot,Y_plot,x_inlnum,y_inlnum,Map_cor,H,pix_tol);
%% plot
if strcmp(dataname, 'boat')
    imfile1 = [dataname '/img' num2str(id1) '.pgm'];
    I1 = imread(imfile1) ;
    imfile2 = [dataname '/img' num2str(id2) '.pgm'];
    I2 = imread(imfile2) ;
    d3 = 1;
else
    imfile1 = [dataname '/img' num2str(id1) '.ppm'];
    I1 = imread(imfile1) ;
    imfile2 = [dataname '/img' num2str(id2) '.ppm'];
    I2 = imread(imfile2) ;
    d3 = 3;
end


size0 = size(I1);
gap = 10;
x_shift = gap + size0(2);

im12 = [I1,uint8(255*ones(size0(1),gap,d3)),I2];
[Dis] = plotMatch_sift(asgHun(Map),XPoint,YPoint,output.X_inl,output.Y_inl,im12,x_shift,H,pix_tol);
