function showPascalGT(fmat, Pts1, Pts2, X_GT)
load(fmat);
im = appendimages(I1, I2);

imshow(im);
hold on;

if nargin < 4
    for i = 1:size(features1, 1)
        r = 5;
        x = features1(i, 2) - r;
        y = features1(i, 1) - r;
        rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
    end

    for i = 1:size(features2, 1)
        r = 5;
        x = features2(i, 2) - r + size(I1, 2);
        y = features2(i, 1) - r;
        rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
    end

    for i = 1:size(gTruth, 2)
        x = [features1(i, 2), features2(gTruth(i), 2) + size(I1, 2)];
        y = [features1(i, 1), features2(gTruth(i), 1)];
        plot(x, y, 'Linewidth', 2, 'Color', 'r');
    end
else
    for i = 1:size(Pts1, 1)
        r = 5;
        x = Pts1(i, 1) - r;
        y = Pts1(i, 2) - r;
        rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
    end

    for i = 1:size(Pts2, 1)
        r = 5;
        x = Pts2(i, 1) - r + size(I1, 2);
        y = Pts2(i, 2) - r;
        rectangle('Position',[x, y, 2 * r, 2 * r], 'Curvature', [1, 1], 'EdgeColor', 'b');
    end
    
    [p, q, ~] = find(X_GT);
    for i = 1:size(p, 1)
        x = [Pts1(p(i), 1), Pts2(q(i), 1) + size(I1, 2)];
        y = [Pts1(p(i), 2), Pts2(q(i), 2)];
        plot(x, y, 'Linewidth', 2, 'Color', 'r');
    end
   
end

hold off;

end