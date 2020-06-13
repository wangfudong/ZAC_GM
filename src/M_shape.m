function [MS,BH1,BH2] = M_shape(X,Y,r0,r1,multi_rota)
LX = length(X(:,1));
LY = length(Y(:,1));

if nargin > 4 && multi_rota > 0
    [BH1] = SC_invariant(X',zeros(1,LX),[],12,5,r0,r1,zeros(1,LX));%BH1(i,:)' is formed form a 5*12 matrix.
    [BH2] = SC_invariant(Y',zeros(1,LY),[],12,5,r0,r1,zeros(1,LY));
    
else
    
    [BH1] = sc_compute(X',zeros(1,LX),[],12,5,r0,r1,zeros(1,LX));%BH1(i,:)' is formed form a 5*12 matrix.
    [BH2] = sc_compute(Y',zeros(1,LY),[],12,5,r0,r1,zeros(1,LY));
end

BH1 = full(BH1);
BH2 = full(BH2);

%     MS = zeros(LX,LY);
%     for i = 1:60
%         MS = MS + abs(bsxfun(@minus,BH1(:,i),BH2(:,i)')).^2;
%     end
%     MS = MS.^(1/2);
%     MS = MS/max(MS(:));

MS = hist_cost_2(BH1,BH2);
MS = MS/max(MS(:));
    
