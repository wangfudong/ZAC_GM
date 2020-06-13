

The ICP algorithm
=================

Per Bergström


This is a short documatation about the ICP (iterative closest point) algorithm implemented in "icp2.m". Examples of its simplest
usage are found in "example.m" and "example2.m" where least-squares point-to-point minimization are used, see [Besl & McKay 1992].
In addition to least-squares minimization other criterion functions are implemented as well. These are:

 1) Huber criterion function (robust)
 2) Tukey's bi-weight criterion function (robust)
 3) Cauchy criterion function (robust)
 4) Welsch criterion function (robust)   

An example where Welsch criterion function is used is found in "example3.m". More information about the robust IRLS-ICP algorithm
is given in [Bergström & Edlund 2014]. See also the documentation about "icp.m" by running "help icp" in Matlab.



Reference:

[Bergström & Edlund 2014]

Bergström, P. and Edlund, O. 2014, “Robust registration of point sets using iteratively reweighted least squares?
Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, doi: 10.1007/s10589-014-9643-2



Further reading:

[Bergström & Edlund (2016) 2017]

Bergström, P. and Edlund, O. (2016) 2017, “Robust registration of surfaces using a refined iterative closest point algorithm 
with a trust region approach? Numerical Algorithms, doi: 10.1007/s11075-016-0170-3


[Bergström 2016]

Bergström, P. 2016, “Reliable updates of the transformation in the iterative closest point algorithm?
Computational Optimization and Applications, vol 63, no. 2, pp. 543-557, doi: doi:10.1007/s10589-015-9771-3


[Bergström, Edlund, & Söderkvist 2011]

Bergström, P. Edlund, O., Söderkvist, I. 2011, “Repeated surface registration for on-line use?
The International Journal of Advanced Manufacturing Technology, vol 54, no. 5-8, pp. 677-689, doi: 10.1007/s00170-010-2950-6



Doi links:

http://dx.doi.org/10.1007/s10589-014-9643-2

http://dx.doi.org/10.1007/s11075-016-0170-3

http://dx.doi.org/10.1007/s10589-015-9771-3

http://dx.doi.org/10.1007/s00170-010-2950-6



Springer Nature’s SharedIt links (full paper online access):

http://rdcu.be/nJRM

http://rdcu.be/noHE

http://rdcu.be/nJRV

http://rdcu.be/nJUW



YouTube: Real time 3D shape analysis by comparing point cloud with CAD-model

Examples where an ICP-algorithm is used can be seen on YouTube,


https://youtu.be/lm7_mwpOk0E


https://youtu.be/cPS-DY9sCz4


A demonstration of shape matching, shape comparison, and shape identification is presented.


The IGES-toolbox is used to read the CAD-models into Matlab,


http://www.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox





