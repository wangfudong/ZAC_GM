 function [dist, X] = FGM_SC(KP, KQ, A1, A2, S1, S2, beta, bDirect)   % fgm
 
    %% FGM_SC cannot work when two few symmetric pairs
    count1 = sum(sum(S1>0));
    count2 = sum(sum(S2>0));
    if count1 < 4 || count2 < 4
        [dist, X] = FGM(KP, KQ, A1, A2, bDirect);
        return;
    end
    
    Cons = ones(size(A1, 1), size(A2, 1));
    
    if ~bDirect
        A1 = tril(A1);
        A2 = tril(A2);
        S1 = tril(S1);
        S2 = tril(S2);
    end
    [p1, q1, v1] = find(A1);
    [p2, q2, v2] = find(A2);
    [sp1, sq1, sv1] = find(S1);
    [sp2, sq2, sv2] = find(S2);
    [n1, n2] = size(KP);
    m1 = size(p1,1);    
    m2 = size(p2,1);
    sm1 = size(sp1,1);
    sm2 = size(sp2,1);
    

    
    %% build matrix used in the FGM algorithm
% %     if ~bDirect
% %         gphs{1}.G = zeros(n1, m1 + sm1);
% %         gphs{2}.G = zeros(n2, m2 + sm2);
% %         for i = 1:m1
% %             gphs{1}.G(p1(i), i) = 1;    
% %             gphs{1}.G(q1(i), i) = 1;
% %         end
% %         for i = 1:sm1
% %             gphs{1}.G(sp1(i), m1 + i) = 1;    
% %             gphs{1}.G(sq1(i), m1 + i) = 1;
% %         end
% %         for i = 1:m2
% %             gphs{2}.G(p2(i), i) = 1;
% %             gphs{2}.G(q2(i), i) = 1;
% %         end
% %         for i = 1:sm2
% %             gphs{2}.G(sp2(i), m2 + i) = 1;
% %             gphs{2}.G(sq2(i), m2 + i) = 1;
% %         end
% %         gphs{1}.H = [gphs{1}.G diag(ones(1, n1))];
% %         gphs{2}.H = [gphs{2}.G diag(ones(1, n2))];
% %     else
% %         gphs{1}.Eg = [p1', sp1'; q1', sq1'];
% %         gphs{2}.Eg = [p2', sp2'; q2', sq2'];
% %         [gphs{1}.G, gphs{1}.H] = gphEg2IncA(gphs{1}.Eg, n1);
% %         [gphs{2}.G, gphs{2}.H] = gphEg2IncA(gphs{2}.Eg, n2);
% %     end
% %     
% %    
% %     KQ_new = zeros(m1 + sm1, m2 + sm2);
% %     for row = 1:m1
% %         for col = 1:m2
% %             KQ_new(row, col) = KQ(row, col);
% %         end
% %     end
% %     
% %     maxv = max([1, max(sv1), max(sv2)]);
% %     for row = 1:sm1
% %         for col = 1:sm2
% %             KQ_new(m1 + row, m2 + col) = beta * (maxv - abs(sv1(row) - sv2(col)));
% %         end
% %     end
% %     
% %     
% %     par = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n');
% %     if bDirect
% %         asg = fgmD(KP, KQ_new, Cons, gphs, [], par);
% %     else
% %         asg = fgmU(KP, KQ_new, Cons, gphs, [], par);
% %     end
% %     X = asg.X;
% %     dist = asg.obj;
    
    
    if ~bDirect
        gphs{1}.G = zeros(n1, m1);
        gphs{2}.G = zeros(n2, m2);
        for i = 1:m1
            gphs{1}.G(p1(i), i) = 1;    
            gphs{1}.G(q1(i), i) = 1;
        end
        for i = 1:m2
            gphs{2}.G(p2(i), i) = 1;
            gphs{2}.G(q2(i), i) = 1;
        end
        gphs{1}.H = [gphs{1}.G diag(ones(1, n1))];
        gphs{2}.H = [gphs{2}.G diag(ones(1, n2))];

        gphs{1}.SG = zeros(n1, sm1);
        gphs{2}.SG = zeros(n2, sm2);
        for i = 1:sm1
            gphs{1}.SG(sp1(i), i) = 1;    
            gphs{1}.SG(sq1(i), i) = 1;
        end
        for i = 1:sm2
            gphs{2}.SG(sp2(i), i) = 1;
            gphs{2}.SG(sq2(i), i) = 1;
        end
        gphs{1}.SH = [gphs{1}.SG diag(ones(1, n1))];
        gphs{2}.SH = [gphs{2}.SG diag(ones(1, n2))];


        KPS = zeros(size(KP));

        KS = zeros(sm1, sm2);
        maxv = max([1, max(sv1), max(sv2)]);
        for row = 1:sm1
            for col = 1:sm2
                KS(row, col) = beta * (maxv - abs(sv1(row) - sv2(col)));
            end
        end

        par = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n');
        asg = fgmU_SC(KP, KQ, KPS, KS, Cons, gphs, [], par);
    else
        gphs{1}.Eg = [p1', sp1'; q1', sq1'];
        gphs{2}.Eg = [p2', sp2'; q2', sq2'];
        [gphs{1}.G, gphs{1}.H] = gphEg2IncA(gphs{1}.Eg, n1);
        [gphs{2}.G, gphs{2}.H] = gphEg2IncA(gphs{2}.Eg, n2);

        KQ_new = zeros(m1 + sm1, m2 + sm2);
        for row = 1:m1
            for col = 1:m2
                KQ_new(row, col) = KQ(row, col);
            end
        end

        maxv = max([1, max(sv1), max(sv2)]);
        for row = 1:sm1
            for col = 1:sm2
                KQ_new(m1 + row, m2 + col) = beta * (maxv - abs(sv1(row) - sv2(col)));
            end
        end


        par = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n');
        asg = fgmD(KP, KQ_new, Cons, gphs, [], par);
   end
    
    X = asg.X;
    dist = asg.obj;
 end
