 function [dist, X] = FGM(KP, KQ, A1, A2, bDirect)   % fgm
    Cons = ones(size(A1, 1), size(A2, 1));
    
    if ~bDirect
        A1 = tril(A1);
        A2 = tril(A2);
    end
    
    [p1,q1,~] = find(A1);
    [p2,q2,~] = find(A2);
    [n1, n2] = size(KP);
    m1 = size(p1,1);    
    m2 = size(p2,1);
    
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
        gphs{1}.H = [gphs{1}.G eye(n1)];
        gphs{2}.H = [gphs{2}.G eye(n2)];
    else
        gphs{1}.Eg = [p1'; q1'];
        gphs{2}.Eg = [p2'; q2'];
        [gphs{1}.G, gphs{1}.H] = gphEg2IncA(gphs{1}.Eg, n1);
        [gphs{2}.G, gphs{2}.H] = gphEg2IncA(gphs{2}.Eg, n2);
    end
   
    
    par = st('nItMa', 100, 'nAlp', 101, 'thAlp', 0, 'deb', 'n');
    if bDirect
        asg = fgmD(KP, KQ, Cons, gphs, [], par);
    else
        asg = fgmU(KP, KQ, Cons, gphs, [], par);
    end
    X = asg.X;
    dist = asg.obj;
 end
