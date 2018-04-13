%the demo of the video background subtraction
%
% References:
% Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali,
% "Low-rank and Sparse Matrix Decomposition via the Truncated Nuclear Norm and a Sparse Regularizer,"
% submitted to The Visual Computer, April 2018.
%
% Written by Zhichao Xue, version 1.0

lrs_setup;

lrs_load_conf;

input_avi = fullfile(lrs_conf.lrs_dir,'dataset','escalator.avi');   %input video, for example: escalator.avi. And Video_003.avi is from BMC dataset.
output_avi = fullfile(lrs_conf.lrs_dir,'output','output.avi');

process_video(input_avi, output_avi);