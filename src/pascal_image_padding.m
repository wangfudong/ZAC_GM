function [X,Y,im12] = pascal_image_padding(im1,im2,X,Y,gap,size0)
size1 = size(im1);
size2 = size(im2);
y1 = floor((size0(1)-size1(1))/2);
y2 = floor((size0(1)-size2(1))/2);
x1 = floor((size0(2)-size1(2))/2);
x2 = floor((size0(2)-size2(2))/2);

im11 = uint8(255*ones(size0));
im22 = uint8(255*ones(size0));

im11(y1+1:y1+size1(1),x1+1:x1+size1(2),:) = im1;
im22(y2+1:y2+size2(1),x2+1:x2+size2(2),:) = im2;

X(:,1) = X(:,1) + x1;
X(:,2) = X(:,2) + y1;
Y(:,1) = Y(:,1) + x2;
Y(:,2) = Y(:,2) + y2;

if length(size1) == 2
    im12 = [im11,uint8(255*ones(size0(1),gap)),im22];
else
    im12 = [im11,uint8(255*ones(size0(1),gap,3)),im22];
end
