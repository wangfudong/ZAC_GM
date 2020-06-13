MATLAB demo code of Max-Pooling Matching CVPR 2014

M. Cho, J. Sun, O. Duchenne, J. Ponce
Finding Matches in a Haystack: A Max-Pooling Strategy for Graph Matching in the Presence of Outliers 
Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (2014) 
http://www.di.ens.fr/willow/research/maxpoolingmatching/

Please cite our work if you find this code useful in your research. 

written by Minsu Cho, Inria - WILLOW / Ecole Normale Superieure 
http://www.di.ens.fr/~mcho/


Date: 20/06/2014
Version: 1.0

==================================================================================================


1. Overview

do_PointMatchingTest.m   : main script for point matching demo

Current test parameters are set to deformation-varing experiments.
For different settings, modify 'setPointMatching.m' in the folder 'utils'.
Refer to our paper for the original setting. 
Performance of each algorithm is represented in terms of accuracy, score,and time.
   Accuracy: the ratio between # of true positive matches and # of groundtruth matches.
   Score: the sum of all affinity values related to the matching result.               


compile.m             : script for c-mex compile

If c-mex functions are incompatible, re-comple c-mex functions by running 'compile.m'

setMethods.m        : script for settings of algorithms being tested
                      the list of methods and their options for visualization

MPM.m              : Matlab function of Max-Pooling Matching 

If you want to add your own algorithm for comparison, three steps are required:
1. Create 'YOUR_ALGORITHM_NAME' folder in 'Methods' folder. Then put your code in it.
2. Add the folder in the script 'setPath.m' so that your method can be called.
3. Modify 'setMethods.m' for your method. Note that you should follow the 'methods' structure. 


2. References

This code includes three algorithms for comparison:

M. Cho, J. Lee, K. M. Lee. "Reweighted Random Walks for Graph Matching", ECCV 2010. 
http://cv.snu.ac.kr/research/~RRWM/
M. Leordeanu, M. Hebert and R. Sutkhankar, "An Integer Projected Fixed Point Method for Graph Matching and MAP Inference", NIPS 2009. 
M. Leordeanu and M. Hebert. "A spectral technique for correspondence problems using pairwise constraints.", ICCV 2005.
https://sites.google.com/site/graphmatchingmethods/

We utilized some functions of the following public implementations;

bistocastic normalization functions of Timothee Cour's: 
http://www.seas.upenn.edu/~timothee/software/graph_matching/graph_matching.html

Hungarian algorithm of Markus Buehren's for final discretization (optional):
http://www.markusbuehren.de/