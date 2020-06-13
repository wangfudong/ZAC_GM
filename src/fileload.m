function [Pts,Ims,nF] = fileload(filenum,inds)
% load the dataset of cmu(filenum=1,2,inds=[g1,g2]) or
% pascal(filenum = 3,4, inds = g12)
if nargin < 2
    error('not enough inputs');
end
if filenum <=2 
    te = length(inds);
    if te ~= 2
        error('the input inds shoule be a 1*2 vector');
    end
else
    te = length(inds);
    if te ~= 1
        error('the input inds should be a integer');
    end
end

Pts = cell(1,2);
Ims = cell(1,2);
if filenum <= 2
    dataset = ['house';'hotel'];
    g1 = inds(1);
    numb1 = num2str(g1,'%03d');
    file1 = [dataset(filenum,:) numb1];
    Pts{1,1} = load(file1);
    Ims{1,1} = imread([dataset(filenum,:) '.seq' num2str(g1-1) '.png']);
    
    g2 = inds(2);
    numb2 = num2str(g2,'%03d');
    file2 = [dataset(filenum,:) numb2];
    Pts{1,2} = load(file2);
    Ims{1,2} = imread([dataset(filenum,:) '.seq' num2str(g2-1) '.png']);
   
    nF = size(Pts{1,1},1);
    
else
    
    file1 = 'Carss/pair_';
    file2 = 'Motor/pair_';
    
    g12 = inds;
    numb1 = num2str(g12);
    if filenum == 3
        file = [file1 numb1];
    else
        file = [file2 numb1];
    end
    pairs = load([file '.mat']);
    
    Pts{1,1} = pairs.features1(:,2:-1:1);
    Pts{1,2} = pairs.features2(:,2:-1:1);
    Ims{1,1} = pairs.I1;
    Ims{1,2} = pairs.I2;
    
    nF = pairs.nF1;
    
end