function [Xnew,Inew] = pascal_graph_padding(X,I,size0)
size1 = size(I);
y1 = floor((size0(1)-size1(1))/2);
x1 = floor((size0(2)-size1(2))/2);

Inew = uint8(255*ones(size0));

Inew(y1+1:y1+size1(1),x1+1:x1+size1(2),:) = I;

Xnew(:,1) = X(:,1) + x1;
Xnew(:,2) = X(:,2) + y1;


