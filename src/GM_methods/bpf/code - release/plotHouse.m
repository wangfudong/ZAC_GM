function plotHouse()

plotResults(20);
%plotResults(30);

end

function plotResults(n1)

[avgAcc, avgObj, avgTm] = calcResults(n1);


fac = 10:10:100;
nNodes = [n1 30];

title = sprintf('%d-%d-Acc',  nNodes(1), nNodes(2));
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
if n1 >= 25
    axis([fac(1), fac(end), 75, 100.1]); 
    set(gca,'Xtick', fac, 'Ytick', [75:5:100], 'Fontsize', 16);
else
    axis([fac(1), fac(end), 25, 85]); 
    set(gca,'Xtick', fac, 'Ytick', [25:10:85], 'Fontsize', 16);
end
legend(['IPFP'], ['PSM'], ['RRWM'],['GAGM'],['FGM'],  ['GNCCP'], ['BPF-G'],...
    'Location', 'Best');
xlabel('gap','FontSize',20)
ylabel('matching accuracy (%)','FontSize',20)
hold off;



maxv = [avgObj.Psm; avgObj.Rrwm; avgObj.Gnccp; avgObj.BpfG; avgObj.FgmU; avgObj.Gagm];
maxv = max(maxv);
title = sprintf('%d-%d-Obj',  nNodes(1), nNodes(2));
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
if n1 >= 30
    axis([fac(1), fac(end), 0.9, 1.001]); 
    set(gca,'Xtick', fac, 'Ytick',  0.9:0.02:1, 'Fontsize', 16);
else
    axis([fac(1), fac(end), 0.8, 1.001]); 
    set(gca,'Xtick', fac, 'Ytick',  0.8:0.04:1, 'Fontsize', 16);
end
legend(['IPFP'], ['PSM'], ['RRWM'], ['GAGM'], ['FGM'],  ['GNCCP'], ['BPF-G'], ...
    'Location', 'Best');
xlabel('gap','FontSize',20)
ylabel('objective ratio','FontSize',20)
hold off;


maxv = ceil(max([max(avgTm.Psm), max(avgTm.FgmU), max(avgTm.Gagm), max(avgTm.Gnccp), max(avgTm.BpfG)]));
minv = 0;
step = (maxv-minv) / 10;
title = sprintf('%d-%d-CpuTime',  nNodes(1), nNodes(2));
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
set(gca,'Xtick', fac, 'Ytick', minv:step:maxv);
legend(['IPFP'], ['PSM'], ['RRWM'],['GAGM'],['FGM'],  ['GNCCP'], ['BPF-G'], ...
    'Location', 'Best');
xlabel('gap','FontSize',16)
ylabel('Computational time','FontSize',16)
hold off;


end


function [avgAcc, avgObj, avgTm] = calcResults(n1)
avgAcc.Psm      = zeros(1, 10);
avgAcc.Gnccp    = zeros(1, 10);
avgAcc.BpfG     = zeros(1, 10);
avgAcc.IpfpS    = zeros(1, 10);
avgAcc.Rrwm     = zeros(1, 10);
avgAcc.FgmU     = zeros(1, 10);
avgAcc.Gagm     = zeros(1, 10);

avgObj = avgAcc;
avgTm  = avgAcc;

savePath = './save/cmum';
for gi = 1:10
    gap = gi * 10;
    for fid = 1:110-gap
        smat_ipfps = sprintf('%s/n%d_g%d_f%d_ipfps.mat', savePath, n1, gap, fid);
        smat_rrwm = sprintf('%s/n%d_g%d_f%d_rrwm.mat', savePath, n1, gap, fid);
        smat_psm = sprintf('%s/n%d_g%d_f%d_psm.mat', savePath, n1, gap, fid);
        smat_gnccp = sprintf('%s/n%d_g%d_f%d_gnccp.mat', savePath, n1, gap, fid);
        smat_bpfG = sprintf('%s/n%d_g%d_f%d_bpfG.mat', savePath, n1, gap, fid);
        smat_fgm = sprintf('%s/n%d_g%d_f%d_fgm.mat', savePath, n1, gap, fid);
        smat_gagm = sprintf('%s/n%d_g%d_f%d_gagm.mat', savePath, n1, gap, fid);
        
        
        load(smat_ipfps);
        avgTm.IpfpS(gi) = avgTm.IpfpS(gi) + asgIpfpS.tm;
        avgAcc.IpfpS(gi) = avgAcc.IpfpS(gi) + asgIpfpS.acc;
        avgObj.IpfpS(gi) = avgObj.IpfpS(gi) + asgIpfpS.obj;
        
        load(smat_rrwm);
        avgTm.Rrwm(gi) = avgTm.Rrwm(gi) + asgRrwm.tm;
        avgAcc.Rrwm(gi) = avgAcc.Rrwm(gi) + asgRrwm.acc;
        avgObj.Rrwm(gi) = avgObj.Rrwm(gi) + asgRrwm.obj;
        
        load(smat_psm);
        avgTm.Psm(gi) = avgTm.Psm(gi) + asgPsm.tm;
        avgAcc.Psm(gi) = avgAcc.Psm(gi) + asgPsm.acc;
        avgObj.Psm(gi) = avgObj.Psm(gi) + asgPsm.obj;
        
        load(smat_gnccp);
        avgTm.Gnccp(gi) = avgTm.Gnccp(gi) + asgGnccp.tm;
        avgAcc.Gnccp(gi) = avgAcc.Gnccp(gi) + asgGnccp.acc;
        avgObj.Gnccp(gi) = avgObj.Gnccp(gi) + asgGnccp.obj;
        
        load(smat_bpfG);
        avgTm.BpfG(gi) = avgTm.BpfG(gi) + asgBpfG.tm;
        avgAcc.BpfG(gi) = avgAcc.BpfG(gi) + asgBpfG.acc;
        avgObj.BpfG(gi) = avgObj.BpfG(gi) + asgBpfG.obj;
        
        load(smat_fgm);
        avgTm.FgmU(gi) = avgTm.FgmU(gi) + asgFgmU.tm;
        avgAcc.FgmU(gi) = avgAcc.FgmU(gi) + asgFgmU.acc;
        avgObj.FgmU(gi) = avgObj.FgmU(gi) + asgFgmU.obj;
        
        load(smat_gagm);
        avgTm.Gagm(gi) = avgTm.Gagm(gi) + asgGagm.tm;
        avgAcc.Gagm(gi) = avgAcc.Gagm(gi) + asgGagm.acc;
        avgObj.Gagm(gi) = avgObj.Gagm(gi) + asgGagm.obj;
        
    end
    
    t = fid;
    
    avgAcc.Psm(gi) = avgAcc.Psm(gi) / t;
    avgAcc.Gnccp(gi) = avgAcc.Gnccp(gi) / t;
    avgAcc.BpfG(gi) = avgAcc.BpfG(gi) / t;
    avgAcc.IpfpS(gi) = avgAcc.IpfpS(gi) / t;
    avgAcc.Rrwm(gi) = avgAcc.Rrwm(gi) / t;
    avgAcc.FgmU(gi) = avgAcc.FgmU(gi) / t;
    avgAcc.Gagm(gi) = avgAcc.Gagm(gi) / t;
    
    avgObj.Psm(gi) = avgObj.Psm(gi) / t;
    avgObj.Gnccp(gi) = avgObj.Gnccp(gi) / t;
    avgObj.BpfG(gi) = avgObj.BpfG(gi) / t;
    avgObj.IpfpS(gi) = avgObj.IpfpS(gi) / t;
    avgObj.Rrwm(gi) = avgObj.Rrwm(gi) / t;
    avgObj.FgmU(gi) = avgObj.FgmU(gi) / t;
    avgObj.Gagm(gi) = avgObj.Gagm(gi) / t;

    avgTm.Psm(gi) = avgTm.Psm(gi) / t;
    avgTm.Gnccp(gi) = avgTm.Gnccp(gi) / t;
    avgTm.BpfG(gi) = avgTm.BpfG(gi) / t;
    avgTm.IpfpS(gi) = avgTm.IpfpS(gi) / t;
    avgTm.Rrwm(gi) = avgTm.Rrwm(gi) / t;
    avgTm.FgmU(gi) = avgTm.FgmU(gi) / t;
    avgTm.Gagm(gi) = avgTm.Gagm(gi) / t;
end

end