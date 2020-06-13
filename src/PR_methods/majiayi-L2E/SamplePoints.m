function [xs1, ys1, ts1, xs2, ys2, ts2] = SamplePoints(S1, nsamp1, S2, nsamp2)
%%%
%    S1 : shape 1
%    nsamp1 : sample count for shape 1
%    S2 : shape 2
%    nsamp2 : sample count for shape 2
%%%

% extract contour
V1 = double(S1);
[x1, y1, t1, c1]=bdry_extract_3(V1);

V2 = double(S2);
[x2, y2, t2, c2]=bdry_extract_3(V2);

% sample points
[xs1, ys1, ts1]=get_samples_1(x1, y1, t1, nsamp1);
[xs2, ys2, ts2]=get_samples_1(x2, y2, t2, nsamp2);

end
