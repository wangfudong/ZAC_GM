function showWillowGT(cdata, X)
im = appendimages(cdata.im1, cdata.im2);
imshow(im);
hold on;

Pts1 = cdata.Pts1';
Pts2 = cdata.Pts2';

%% draw points
for i = 1:size(Pts1, 1)
    r = 5;
    x = Pts1(i, 1) - r;
    y = Pts1(i, 2) - r;
    rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
end

width = size(cdata.im1, 2);

for i = 1:size(Pts2, 1)
    r = 5;
    x = Pts2(i, 1) - r + width;
    y = Pts2(i, 2) - r;
    rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
end

%% draw edges
[p,q,~] = find(tril(cdata.A1));
for i = 1:size(p,1)
    x = [Pts1(p(i),1), Pts1(q(i),1)];
    y = [Pts1(p(i),2), Pts1(q(i),2)];
    plot(x, y, 'LineWidth', 2, 'Color', 'y');
end


[p,q,~] = find(tril(cdata.A2));
for i = 1:size(p,1)
    x = [Pts2(p(i),1), Pts2(q(i),1)] + width;
    y = [Pts2(p(i),2), Pts2(q(i),2)];
    plot(x, y, 'LineWidth', 2, 'Color', 'y');
end


%% draw ground truth
if size(X, 2) > 1
    [p, q, ~] = find(X);
    for i = 1:size(p, 1)
        x = [Pts1(p(i), 1), Pts2(q(i), 1) + width];
        y = [Pts1(p(i), 2), Pts2(q(i), 2)];
        
        if cdata.X_GT(p(i), q(i))
            clr = 'g';
        else
            clr = 'r';
        end
        
        plot(x, y, 'Linewidth', 2, 'Color', clr);
    end
else
    ind = find(X);
    for i = 1:size(ind, 1)
        p = cdata.matchInfo(1, ind(i));
        q = cdata.matchInfo(2, ind(i));
        x = [Pts1(p, 1), Pts2(q, 1) + width];
        y = [Pts1(p, 2), Pts2(q, 2)];
        
        if cdata.x_gt(ind(i))
            clr = 'g';
        else
            clr = 'r';
        end
        
        plot(x, y, 'Linewidth', 2, 'Color', clr);
    end        
end


hold off;

end