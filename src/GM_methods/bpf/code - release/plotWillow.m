function plotWillow()

savePath = './save/Willow';
[avgAcc, avgObj, avgCT] = getResult(savePath);

fac = 0:10;

title = 'Willow-Acc';
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
legend( ['IPFP'], ['PSM'], ['RRWM'], ['GAGM'], ['FGM'],['GNCCP'], ['BPF-G'],...
     'Location', 'Best');
axis([fac(1), fac(end), 30, 100]); 
set(gca,'Xtick', fac, 'Ytick', [30:10:100], 'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('matching accuracy (%)','FontSize',20)
hold off;



maxv = [avgObj.Psm; avgObj.Rrwm; avgObj.Gnccp; avgObj.BpfG; avgObj.FgmU];
maxv = max(maxv);
title = 'Willow-Obj';
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
legend(['IPFP'], ['PSM'], ['RRWM'], ['GAGM'],['FGM'], ['GNCCP'], ['BPF-G'], ...
     'Location', 'Best');
axis([fac(1), fac(end), 0.85, 1.001]); 
set(gca,'Xtick', fac, 'Ytick', 0.85:0.05:1, 'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('objective ratio','FontSize',20)
hold off;


maxv = ceil(max([max(avgCT.Psm), max(avgCT.Gnccp), max(avgCT.BpfG), max(avgCT.FgmU)]));
minv = 0;
step = (maxv-minv) / 10;
title = 'Willow-CpuTime';
figure('NumberTitle', 'off', 'Name', title);
hold on;
plot(fac, avgCT.IpfpS, 'k*--', ...
     fac, avgCT.Psm, 'co--', ...
     fac, avgCT.Rrwm, 'gv--', ...
     fac, avgCT.Gagm, 'bx--', ...
     fac, avgCT.FgmU, 'k+-', ...
     fac, avgCT.Gnccp, 'rd-', ...
     fac, avgCT.BpfG, 'ms-', ...
     'LineWidth', 2, 'MarkerSize', 10); 
legend(['IPFP'], ['PSM'], ['RRWM'],  ['GAGM'],['FGM'],['GNCCP'], ['BPF-G'],  ...
     'Location', 'Best');
axis([fac(1), fac(end), minv, maxv]); 
set(gca,'Xtick', fac, 'Ytick', minv:step:maxv, 'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('Computational time','FontSize',20)
hold off;
end


function [avgAcc, avgObj, avgCT] = getResult(savePath)
outliers = 0:2:20;

avgAcc.Psm      = zeros(size(outliers));
avgAcc.Gnccp    = zeros(size(outliers));
avgAcc.BpfG     = zeros(size(outliers));
avgAcc.IpfpS    = zeros(size(outliers));
avgAcc.Rrwm     = zeros(size(outliers));
avgAcc.FgmU     = zeros(size(outliers));
avgAcc.Gagm     = zeros(size(outliers));

avgObj = avgAcc;
avgCT  = avgAcc;

imgClass = {'Car', 'Duck', 'Face', 'Motorbike', 'Winebottle'};
t = 0;
for cls = 1:size(imgClass, 2)
    clsName = imgClass{cls};
    
    for cImg = 1:100
        t = t + 1;
        
        smat_ipfps  = sprintf('%s/%s_%04d_ipfpS.mat', savePath, clsName, cImg);
        smat_psm    = sprintf('%s/%s_%04d_psm.mat', savePath, clsName, cImg);
        smat_rrwm   = sprintf('%s/%s_%04d_rrwm.mat', savePath, clsName, cImg);
        smat_gnccp  = sprintf('%s/%s_%04d_gnccp.mat', savePath, clsName, cImg);
        smat_BpfG = sprintf('%s/%s_%04d_bpfg.mat', savePath, clsName, cImg);
        smat_fgm    = sprintf('%s/%s_%04d_fgm.mat', savePath, clsName, cImg);
        smat_gagm    = sprintf('%s/%s_%04d_gagm.mat', savePath, clsName, cImg);
        
        load(smat_ipfps);
        load(smat_psm);
        load(smat_rrwm);
        load(smat_gnccp);
        load(smat_BpfG);
        load(smat_fgm);
        load(smat_gagm);

       for i = 1:size(outliers,2)
           avgCT.Psm(i)  = avgCT.Psm(i)  + asgPsm{i}.tm;
           avgAcc.Psm(i) = avgAcc.Psm(i) + asgPsm{i}.acc;
           avgObj.Psm(i) = avgObj.Psm(i) + asgPsm{i}.obj;

           avgCT.Gnccp(i) = avgCT.Gnccp(i) + asgGnccp{i}.tm;
           avgAcc.Gnccp(i) = avgAcc.Gnccp(i) + asgGnccp{i}.acc;
           avgObj.Gnccp(i) = avgObj.Gnccp(i) + asgGnccp{i}.obj;

           avgCT.BpfG(i) = avgCT.BpfG(i) + asgBpfG{i}.tm;
           avgAcc.BpfG(i) = avgAcc.BpfG(i) + asgBpfG{i}.acc;
           avgObj.BpfG(i) = avgObj.BpfG(i) + asgBpfG{i}.obj;

           avgCT.IpfpS(i) = avgCT.IpfpS(i) + asgIpfpS{i}.tm;
           avgAcc.IpfpS(i) = avgAcc.IpfpS(i) + asgIpfpS{i}.acc;
           avgObj.IpfpS(i) = avgObj.IpfpS(i) + asgIpfpS{i}.obj;

           avgCT.Rrwm(i) = avgCT.Rrwm(i) + asgRrwm{i}.tm;
           avgAcc.Rrwm(i) = avgAcc.Rrwm(i) + asgRrwm{i}.acc;
           avgObj.Rrwm(i) = avgObj.Rrwm(i) + asgRrwm{i}.obj;

           avgCT.FgmU(i) = avgCT.FgmU(i) + asgFgmD{i}.tm;
           avgAcc.FgmU(i) = avgAcc.FgmU(i) + asgFgmD{i}.acc;
           avgObj.FgmU(i) = avgObj.FgmU(i) + asgFgmD{i}.obj;

           avgCT.Gagm(i) = avgCT.Gagm(i) + asgGagm{i}.tm;
           avgAcc.Gagm(i) = avgAcc.Gagm(i) + asgGagm{i}.acc;
           avgObj.Gagm(i) = avgObj.Gagm(i) + asgGagm{i}.obj;
       end  
    end
end

avgAcc.Psm = avgAcc.Psm / t;
avgAcc.Gnccp = avgAcc.Gnccp / t;
avgAcc.BpfG =  avgAcc.BpfG / t;
avgAcc.IpfpS = avgAcc.IpfpS / t;
avgAcc.Rrwm = avgAcc.Rrwm / t;
avgAcc.FgmU = avgAcc.FgmU / t;
avgAcc.Gagm = avgAcc.Gagm / t;

avgObj.Psm = avgObj.Psm / t;
avgObj.Gnccp = avgObj.Gnccp / t;
avgObj.BpfG = avgObj.BpfG / t;
avgObj.IpfpS = avgObj.IpfpS / t;
avgObj.Rrwm = avgObj.Rrwm / t;
avgObj.FgmU = avgObj.FgmU / t;
avgObj.Gagm = avgObj.Gagm / t;

avgCT.Psm = avgCT.Psm / t;
avgCT.Gnccp = avgCT.Gnccp / t;
avgCT.BpfG = avgCT.BpfG / t;
avgCT.IpfpS = avgCT.IpfpS / t;
avgCT.Rrwm = avgCT.Rrwm / t;
avgCT.FgmU = avgCT.FgmU / t;
avgCT.Gagm = avgCT.Gagm / t;

end