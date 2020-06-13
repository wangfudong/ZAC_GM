function plotPascal()
savePath = './save/Pascal';
[avgAcc, avgObj, avgCT] = getResult(savePath);

fac = 0:2:20;

title = 'Pascal-Acc';
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
legend(['IPFP'], ['PSM'], ['RRWM'], ['GAGM'],['FGM'], ['GNCCP'], ['BPF-G'], ...
     'Location', 'Best');
axis([fac(1), fac(end), 10, 90]); 
set(gca,'Xtick', fac, 'Ytick', [10:10:90], 'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('matching accuracy (%)','FontSize',20)
hold off;



maxv = [avgObj.Psm; avgObj.Rrwm; avgObj.Gnccp; avgObj.BpfG; avgObj.FgmU; avgObj.Gagm]; 
maxv = max(maxv);
title = 'Pascal-Obj';
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
legend(['IPFP'], ['PSM'], ['RRWM'], ['GAGM'], ['FGM'],['GNCCP'], ['BPF-G'],  ...
     'Location', 'Best');
axis([fac(1), fac(end), 0.75, 1.001]); 
set(gca,'Xtick', fac, 'Ytick', 0.75:0.05:1, 'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('objective ratio','FontSize',20)
hold off;


maxv = ceil(max([max(avgCT.Psm), max(avgCT.FgmU), max(avgCT.Gagm), max(avgCT.Gnccp), max(avgCT.BpfG)]));
minv = 0;
step = (maxv-minv) / 10;
title = 'Pascal-CpuTime';
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
legend(['IPFP'], ['PSM'], ['RRWM'],['GAGM'],['FGM'], ['GNCCP'], ['BPF-G'],  ...
     'Location', 'Best');
axis([fac(1), fac(end), minv, maxv]); 
set(gca,'Xtick', fac, 'Ytick', minv:step:maxv,'FontSize', 16);
xlabel('outliers','FontSize',20)
ylabel('Computational time','FontSize',20)
hold off;

end


function [avgAcc, avgObj, avgCT] = getResult(savePath)
avgAcc.Psm      = zeros(1, 11);
avgAcc.Gnccp    = zeros(1, 11);
avgAcc.BpfG     = zeros(1, 11);
avgAcc.IpfpS    = zeros(1, 11);
avgAcc.Rrwm     = zeros(1, 11);
avgAcc.FgmU     = zeros(1, 11);
avgAcc.Gagm     = zeros(1, 11);

avgObj = avgAcc;
avgCT  = avgAcc;

t = 0;
for cImg = 1:50
    t = t + 1;
    if cImg <= 20
        smat_ipfps = sprintf('%s/Motor%02d_ipfps.mat', savePath, cImg);
        smat_psm = sprintf('%s/Motor%02d_psm.mat', savePath, cImg);
        smat_rrwm = sprintf('%s/Motor%02d_rrwm.mat', savePath, cImg);
        smat_gnccp = sprintf('%s/Motor%02d_gnccp.mat', savePath, cImg);
        smat_BpfG = sprintf('%s/Motor%02d_bpfG.mat', savePath, cImg);
        smat_fgm = sprintf('%s/Motor%02d_fgm.mat', savePath, cImg);
        smat_gagm = sprintf('%s/Motor%02d_gagm.mat', savePath, cImg);
    else
        smat_ipfps = sprintf('%s/Car%02d_ipfps.mat', savePath, cImg - 20);
        smat_psm = sprintf('%s/Car%02d_psm.mat', savePath, cImg - 20);
        smat_rrwm = sprintf('%s/Car%02d_rrwm.mat', savePath, cImg - 20);
        smat_gnccp = sprintf('%s/Car%02d_gnccp.mat', savePath, cImg - 20);
        smat_BpfG = sprintf('%s/Car%02d_bpfG.mat', savePath, cImg - 20);
        smat_fgm = sprintf('%s/Car%02d_fgm.mat', savePath, cImg - 20);
        smat_gagm = sprintf('%s/Car%02d_gagm.mat', savePath, cImg - 20);
    end    
    load(smat_ipfps);
    load(smat_psm);
    load(smat_rrwm);
    load(smat_gnccp);
    load(smat_BpfG);
    load(smat_fgm);
    load(smat_gagm);
    
   for i = 1:11
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
       
       avgCT.FgmU(i) = avgCT.FgmU(i) + asgFgm{i}.tm;
       avgAcc.FgmU(i) = avgAcc.FgmU(i) + asgFgm{i}.acc;
       avgObj.FgmU(i) = avgObj.FgmU(i) + asgFgm{i}.obj;

       avgCT.Gagm(i) = avgCT.Gagm(i) + asgGagm{i}.tm;
       avgAcc.Gagm(i) = avgAcc.Gagm(i) + asgGagm{i}.acc;
       avgObj.Gagm(i) = avgObj.Gagm(i) + asgGagm{i}.obj;       
   end  
end

avgAcc.Psm = avgAcc.Psm / t;
avgAcc.Gnccp = avgAcc.Gnccp / t;
avgAcc.BpfG = avgAcc.BpfG / t;
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