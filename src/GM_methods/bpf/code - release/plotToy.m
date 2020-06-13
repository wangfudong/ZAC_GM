function plotToy()

%plotResults('outliers');
%plotResults('noise');
plotResults('density');

end

function plotResults(setting)

if strcmp(setting, 'outliers')
    fac = 0:1:10;
    acc = 0:5:50;
elseif strcmp(setting, 'noise')
    fac = 0:0.02:0.2;
    acc = 0:10:60;
elseif strcmp(setting, 'density')
    fac = 0.3:0.1:1;
    acc = 0:10:100;
end

[avgAcc, avgObj, avgTm] = calcResults(setting);


title = sprintf('%s-Acc',  setting);
figure('NumberTitle', 'off', 'Name', title);
hold on;
plot(fac, 100 * avgAcc.IpfpS, 'k*--', ...
     fac, 100 * avgAcc.Psm, 'co--', ...
     fac, 100 * avgAcc.Rrwm, 'gv--', ...
     fac, 100 * avgAcc.Gagm, 'bx--', ...
     fac, 100 * avgAcc.FgmU, 'k+-', ...
     fac, 100 * avgAcc.Gnccp, 'rd-', ...
     fac, 100 * avgAcc.BpfG, 'ms-', ...
     'LineWidth', 2, 'MarkerSize', 10); 
axis([fac(1), fac(end), acc(1), acc(end) + 1]); 
set(gca,'Xtick', fac, 'Ytick', acc, 'FontSize', 16);
legend(['IPFP'], ['PSM'], ['RRWM'],['GAGM'],['FGM'],['GNCCP'], ['BPF-G'], ...
    'Location', 'Best');
xlabel(setting,'FontSize',20)
ylabel('matching accuracy (%)','FontSize',20)
hold off;



maxv = [avgObj.Psm; avgObj.Rrwm; avgObj.Gnccp; avgObj.BpfG; avgObj.FgmU; avgObj.Gagm];
maxv = max(maxv);
title = sprintf('%s-Obj',  setting);
figure('NumberTitle', 'off', 'Name', title);
hold on;
plot(fac, avgObj.IpfpS ./ maxv, 'k*--', ...
     fac, avgObj.Psm ./ maxv, 'co--', ...
     fac, avgObj.Rrwm ./ maxv, 'gv--', ...
     fac, avgObj.Gagm ./ maxv, 'bx--', ...
     fac, avgObj.FgmU ./ maxv, 'k+-', ...
     fac, avgObj.Gnccp ./ maxv, 'rd-', ...
     fac, avgObj.BpfG ./ maxv, 'ms-', ...
     'LineWidth', 2, 'MarkerSize', 10); 
axis([fac(1), fac(end), 0.6, 1.001]); 
set(gca,'Xtick', fac, 'Ytick',  0.6:0.1:1 , 'FontSize', 16);
legend(['IPFP'], ['PSM'],['RRWM'], ['GAGM'],['FGM'], ['GNCCP'],['BPF-G'], ...
    'Location', 'Best');
xlabel(setting,'FontSize',20)
ylabel('objective ratio','FontSize',20)
hold off;


maxv = ceil(max([max(avgTm.Psm), max(avgTm.FgmU), max(avgTm.Gagm), max(avgTm.BpfG), max(avgTm.Gnccp)]));
minv = 0;
step = (maxv-minv) / 10;
title = sprintf('%s-CpuTime',  setting);
figure('NumberTitle', 'off', 'Name', title);
hold on;
plot(fac, avgTm.IpfpS, 'k*--', ...
     fac, avgTm.Psm, 'co--', ...
     fac, avgTm.Rrwm, 'gv--', ...
     fac, avgTm.Gagm, 'bx--', ...
     fac, avgTm.FgmU, 'k+-', ...
     fac, avgTm.Gnccp, 'rd-', ...
     fac, avgTm.BpfG, 'ms-', ...
     'LineWidth', 2, 'MarkerSize', 10); 
axis([fac(1), fac(end), minv, maxv]);  
set(gca,'Xtick', fac, 'Ytick', minv:step:maxv, 'FontSize', 16);
legend(['IPFP'], ['PSM'], ['RRWM'], ['GAGM'],['FGM'], ['GNCCP'], ['BPF-G'], ...
    'Location', 'Best');
xlabel(setting,'FontSize',20)
ylabel('Computational time','FontSize',20)
hold off;


end


function [avgAcc, avgObj, avgTm] = calcResults(setting)
if strcmp(setting, 'outliers')
    var = 0:1:10;
elseif strcmp(setting, 'noise')
    var = 0:2:20;
elseif strcmp(setting, 'density')
    var = 3:1:10;
end


avgAcc.Psm      = zeros(1, size(var,2));
avgAcc.Gnccp    = zeros(1, size(var,2));
avgAcc.BpfG     = zeros(1, size(var,2));
avgAcc.IpfpS    = zeros(1, size(var,2));
avgAcc.Rrwm     = zeros(1, size(var,2));
avgAcc.FgmU     = zeros(1, size(var,2));
avgAcc.Gagm     = zeros(1, size(var,2));

avgObj = avgAcc;
avgTm  = avgAcc;

savePath = ['./save/toy/' setting];
 
for index = 1:size(var,2)
    v = var(index);
    for gid = 1:100
        smat_fgm  = sprintf('%s/%s%02d_%d_fgm.mat', savePath,setting, v, gid);
        smat_ipfps  = sprintf('%s/%s%02d_%d_ipfps.mat', savePath,setting, v, gid);
        smat_rrwm   = sprintf('%s/%s%02d_%d_rrwm.mat', savePath,setting, v, gid);
        smat_psm    = sprintf('%s/%s%02d_%d_psm.mat', savePath, setting,v, gid);
        smat_gnccp  = sprintf('%s/%s%02d_%d_gnccp.mat', savePath,setting, v, gid);
        smat_BpfG = sprintf('%s/%s%02d_%d_BpfG.mat', savePath, setting,v, gid);
        smat_gagm = sprintf('%s/%s%02d_%d_gagm.mat', savePath, setting,v, gid);

        
        load(smat_fgm);
        avgTm.FgmU(index) = avgTm.FgmU(index) + asgFgmD.tm;
        avgAcc.FgmU(index) = avgAcc.FgmU(index) + asgFgmD.acc;
        avgObj.FgmU(index) = avgObj.FgmU(index) + asgFgmD.obj;
        
        load(smat_ipfps);
        avgTm.IpfpS(index) = avgTm.IpfpS(index) + asgIpfpS.tm;
        avgAcc.IpfpS(index) = avgAcc.IpfpS(index) + asgIpfpS.acc;
        avgObj.IpfpS(index) = avgObj.IpfpS(index) + asgIpfpS.obj;
        
        load(smat_rrwm);
        avgTm.Rrwm(index) = avgTm.Rrwm(index) + asgRrwm.tm;
        avgAcc.Rrwm(index) = avgAcc.Rrwm(index) + asgRrwm.acc;
        avgObj.Rrwm(index) = avgObj.Rrwm(index) + asgRrwm.obj;
        
        load(smat_psm);
        avgTm.Psm(index) = avgTm.Psm(index) + asgPsm.tm;
        avgAcc.Psm(index) = avgAcc.Psm(index) + asgPsm.acc;
        avgObj.Psm(index) = avgObj.Psm(index) + asgPsm.obj;

        load(smat_gnccp);
        avgTm.Gnccp(index) = avgTm.Gnccp(index) + asgGnccp.tm;
        avgAcc.Gnccp(index) = avgAcc.Gnccp(index) + asgGnccp.acc;
        avgObj.Gnccp(index) = avgObj.Gnccp(index) + asgGnccp.obj;
        
        load(smat_BpfG);
        avgTm.BpfG(index) = avgTm.BpfG(index) + asgBpfG.tm;
        avgAcc.BpfG(index) = avgAcc.BpfG(index) + asgBpfG.acc;
        avgObj.BpfG(index) = avgObj.BpfG(index) + asgBpfG.obj;
        
        load(smat_gagm);
        avgTm.Gagm(index) = avgTm.Gagm(index) + asgGagm.tm;
        avgAcc.Gagm(index) = avgAcc.Gagm(index) + asgGagm.acc;
        avgObj.Gagm(index) = avgObj.Gagm(index) + asgGagm.obj;
        
    end
    
    t = gid;
    
    avgAcc.FgmU(index) = avgAcc.FgmU(index) / t;
    avgAcc.Psm(index) = avgAcc.Psm(index) / t;
    avgAcc.Gnccp(index) = avgAcc.Gnccp(index) / t;
    avgAcc.BpfG(index) = avgAcc.BpfG(index) / t;
    avgAcc.IpfpS(index) = avgAcc.IpfpS(index) / t;
    avgAcc.Rrwm(index) = avgAcc.Rrwm(index) / t;
    avgAcc.Gagm(index) = avgAcc.Gagm(index) / t;
   
    avgObj.FgmU(index) = avgObj.FgmU(index) / t;
    avgObj.Psm(index) = avgObj.Psm(index) / t;
    avgObj.Gnccp(index) = avgObj.Gnccp(index) / t;
    avgObj.BpfG(index) = avgObj.BpfG(index) / t;
    avgObj.IpfpS(index) = avgObj.IpfpS(index) / t;
    avgObj.Rrwm(index) = avgObj.Rrwm(index) / t;
    avgObj.Gagm(index) = avgObj.Gagm(index) / t;

    avgTm.FgmU(index) = avgTm.FgmU(index) / t;
    avgTm.Psm(index) = avgTm.Psm(index) / t;
    avgTm.Gnccp(index) = avgTm.Gnccp(index) / t;
    avgTm.BpfG(index) = avgTm.BpfG(index) / t;
    avgTm.IpfpS(index) = avgTm.IpfpS(index) / t;
    avgTm.Rrwm(index) = avgTm.Rrwm(index) / t;
    avgTm.Gagm(index) = avgTm.Gagm(index) / t;
end

end