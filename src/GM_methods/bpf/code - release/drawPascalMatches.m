function drawPascalMatches(fmat, gphs1, gphs2, X, X_GT)
load(fmat);
im = appendimages(I1, I2, 'h', 5);

imshow(im);
hold on;

xshift = 5 + size(I1, 2);
%% draw graph edges
for i = 1:size(gphs1.Eg, 2)
    p1 = gphs1.Eg(1, i);
    p2 = gphs1.Eg(2, i);
    x = [gphs1.Pts(p1, 1), gphs1.Pts(p2, 1)];
    y = [gphs1.Pts(p1, 2), gphs1.Pts(p2, 2)];
    plot(x, y, 'Linewidth', 2, 'Color', 'y');
end

for i = 1:size(gphs2.Eg, 2)
    p1 = gphs2.Eg(1, i);
    p2 = gphs2.Eg(2, i);
    x = [gphs2.Pts(p1, 1) + xshift, gphs2.Pts(p2, 1) + xshift];
    y = [gphs2.Pts(p1, 2), gphs2.Pts(p2, 2)];
    plot(x, y, 'Linewidth', 2, 'Color', 'y');
end
  
%% draw matches
[p, q, ~] = find(X);
for i = 1:size(p, 1)
    if X_GT(p(i), q(i))
        clr = 'g';
    else
        clr = 'r';
    end
    x = [gphs1.Pts(p(i), 1), gphs2.Pts(q(i), 1) + xshift];
    y = [gphs1.Pts(p(i), 2), gphs2.Pts(q(i), 2)];
    plot(x, y, 'Linewidth', 2, 'Color', clr);
end

hold off;

end