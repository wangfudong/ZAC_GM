Introduction
============

This page contains software and instructions for [Branching Path Following for Graph Matching (BPF)] [1] [2].  In addition, we include the following
state-of-the-arts methods as baselines:

- [graduated assignment (GA)] [3]
- [integer projected fixed point method (IPFP)] [4] 
- [re-weighted random walk matching (RRWM)] [5].
- [probabilistic spectral graph matching (PSM)] [6].
- [factorized graph matching (FGM)] [7].
- [graduated non-convexity and concavity procedure(GNCCP)] [8].

The implementations of the above methods are taken from the authors'
websites. We appreciate all the authors for their generosity in sharing codes.


Installation
============

1. unzip `bpf.rar` to your folder;
2. unzip `data.rar` to the subfolder;
3. Run `make` to compile all C++ files;
4. Run `addPath` to add sub-directories into the path of Matlab.
5. Run `demoXXX` or `testXXX`.
6. Run `plotXXX` to show statistical results.


References
==========

[1] T. Wang, H. Ling, C. Lang, and J. Wu. Branching path following for graph matching. In Proceedings of 14th European Conference on Computer Vision (ECCV), 2016. 

[2] T. Wang, H. Ling, C. Lang, and S. Feng. Branching and adaptive path following for graph matching. IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI). (Under Review) 

[3] S. Gold and A. Rangarajan, "A Graduated Assignment Algorithm for Graph Matching", IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 1996.

[4] M. Leordeanu, M. Hebert and R. Sukthankar, "An Integer Projected Fixed Point Method for Graph Matching and MAP Inference", in Advances in Neural Information Processing Systems (NIPS), 2009.

[5] M. Cho, J. Lee and K. Lee, "Reweighted Random Walks for Graph Matching", in European Conference on Computer Vision (ECCV), 2010.

[6] A. Egozi, Y. Keller, and H. Guterman. A probabilistic approach to spectral graph matching. IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 2013.

[7] F. Zhou and F. De la Torre, "Factorized Graph Matching," in IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2012.

[8] Z. Liu and H. Qiao. GNCCP - graduated nonconvexityand concavity procedure.  IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 2014.


Copyright
=========

This software is free for use in research projects. If you publish results obtained using this software, please use this citation.

@inproceedings{WangLLW16,
  author    = {Tao Wang and   Haibin Ling and  Congyan Lang and  Jun Wu},
  title     = {Branching Path Following for Graph Matching},
  booktitle = {European Conference on Computer Vision  (ECCV)},
  pages     = {508--523},
  year      = {2016},
}

If you have any question, please feel free to contact Tao Wang (twang@bjtu.edu.cn).
