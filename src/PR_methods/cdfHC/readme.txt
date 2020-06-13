Group-wise Point-set registration using a novel CDF-based Havrda-Charvat Divergence (CDF-HC) 2D_Demo (version 20111209):
-------------------------------------------------------------------
Copyright (C) 2006-2011 Ting Chen, Baba C. Vemuri and Anand Rangarajan

Authors: Ting Chen
Date:    12/09/2011
 
Contact Information: Ting Chen, tichen@cise.ufl.edu
 
Terms: The source code (M-files) are provided under the
    terms of the GNU General Public License with an explicit
    clause permitting the execution of the M-files from within
    a MATLAB environment. See the LICENSE file for details.
------------------------------------------------------------------

Please cite this paper if you use any pieces of the code.

Reference: Ting Chen, Baba C. Vemuri, Anand Rangarajan and Stephan J. Eisenschenk, 
            Group-wise Point-set registration using a novel CDF-based Havrda-Charvat Divergence. 
            In IJCV : International Journal of Computer Vision, 86(1):111-124, January, 2010.  

--------------------------------------------------------------------

DETAILS:

This directory contains the MATLAB code for the 2d version of the algorithm of group wise point set registration using Havrda-Charvat divergence discribed in the IJCV'10 paper:

"Ting Chen, Baba C. Vemuri, Anand Rangarajan and Stephan J. Eisenschenk, Group-wise Point-set registration using a novel CDF-based Havrda-Charvat Divergence. In IJCV : International Journal of Computer Vision, Volume 86, Number 1, Page 111-124, January, 2010."

Files in this directory are organized as follows: 

	Demo/
		HC2Reg_TPS.m is the main function to start with.
	Data/
		It contains the data sets for testing the program.
	Algorithms/
		This folder contains the implementations for CDF-based Havrda-Charvat Divergence algorithm and some supporting functions for thin plate spline, transformation and display. 
		
Use 'addpath(genpath(pwd))' to add the whole directory and its subdirectories to the MATLAB search path.

Please addresse any questions to tichen@cise.ufl.edu

-----------------------------------------------------
EXAMPLES:

clear all;
close all;
addpath(genpath(pwd));

load cdfhc_data2D_CC
HC2Reg_TPS(CC7,1,300)

%load cdfhc_data2D_beijing
%HC2Reg_TPS(beijing,1)

%load cpd_data2D_fish
%HC2Reg_TPS(cpdfish,1,250)


