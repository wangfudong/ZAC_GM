cd utils
mex assignmentoptimal.cpp
cd ..

% RRWM
cd Methods
cd RRWM
mex mexBistocNormalize_match_slack.cpp 
cd ..
cd MPM
mex -largeArrayDims RMP_mult.cpp 
cd ..
cd ..