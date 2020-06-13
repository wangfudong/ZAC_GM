function drawHouseMatches(t1, t2, gphs1, gphs2, X, X_GT)
imf1 = sprintf('./data/cmum/house/images/house.seq%d.png', t1);
imf2 = sprintf('./data/cmum/house/images/house.seq%d.png', t2);

I1 = imread(imf1);
I2 = imread(imf2);
im = appendimages(I1, I2, 'h', 5);

imshow(im);
hold on;

xshift = 5 + size(I1, 2);
%% draw graph edges
for i = 1:size(gphs1.Eg, 2)
    p1 = gphs1.Eg(1, i);
    p2 = gphs1.Eg(2, i);
    x = [gphs1.Pt(1, p1), gphs1.Pt(1, p2)];
    y = [gphs1.Pt(2, p1), gphs1.Pt(2, p2)];
    plot(x, y, 'Linewidth', 2, 'Color', 'y');
end

for i = 1:size(gphs2.Eg, 2)
    p1 = gphs2.Eg(1, i);
    p2 = gphs2.Eg(2, i);
    x = [gphs2.Pt(1, p1) + xshift, gphs2.Pt(1, p2) + xshift];
    y = [gphs2.Pt(2, p1), gphs2.Pt(2, p2)];
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
    x = [gphs1.Pt(1, p(i)), gphs2.Pt(1, q(i)) + xshift];
    y = [gphs1.Pt(2, p(i)), gphs2.Pt(2, q(i))];
    plot(x, y, 'Linewidth', 2, 'Color', clr);
end


hold off;

end