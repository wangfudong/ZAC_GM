function drawWillowMatches(I1, I2, gphs1, gphs2, X, X_GT)
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
  
%% draw false matches
[p, q, ~] = find(X);
for i = 1:size(p, 1)
    if ~X_GT(p(i), q(i))
        clr = 'r';
        width = 2;
        x = [gphs1.Pts(p(i), 1), gphs2.Pts(q(i), 1) + xshift];
        y = [gphs1.Pts(p(i), 2), gphs2.Pts(q(i), 2)];
        plot(x, y, 'Linewidth', width, 'Color', clr);
    end

end

for i = 1:size(p, 1)
    if X_GT(p(i), q(i))
        clr = 'g';
        width = 3;

        x = [gphs1.Pts(p(i), 1), gphs2.Pts(q(i), 1) + xshift];
        y = [gphs1.Pts(p(i), 2), gphs2.Pts(q(i), 2)];
        plot(x, y, 'Linewidth', width, 'Color', clr);
    end
end


% Pts1 = Pts1';
% Pts2 = Pts2';
% Pts2(:,1) = Pts2(:,1) + size(I1,2) + 5;
% 
% [p1, p2] = find(X);
% 
% % draw matches
% for i = 1:size(p1, 1)
%     if X_GT(p1(i), p2(i)) == 0
%         col1 = 'r';
%         col2 = 'r';
%         linewidth = 2;
%     else
%         col1 = 'g';
%         col2 = 'g';
%         linewidth = 3;
%     end
%     
%     
%     plot([ Pts1(p1(i),1), Pts2(p2(i),1) ]...
%     ,[ Pts1(p1(i),2), Pts2(p2(i),2) ],...
%             '-','LineWidth',linewidth,'MarkerSize',10,...
%             'color', col1);
% 
%     radius = 5;
%     
%     width = 2 * radius;
%     x1 = Pts1(p1(i),1) - radius;
%     y1 = Pts1(p1(i),2) - radius;
%     rectangle('Curvature', [1 1], 'Position', [x1, y1, width, width], 'EdgeColor', col2, 'LineWidth', linewidth);
% 
%     width = 2 * radius;
%     x2 = Pts2(p2(i),1) - radius;
%     y2 = Pts2(p2(i),2) - radius;
%     rectangle('Curvature', [1 1], 'Position', [x2, y2, width, width], 'EdgeColor', col2, 'LineWidth', linewidth);
%     
% end

hold off;

end