%the demo of the face image shadow removal
%
% References:
% Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali,
% "Low-rank and Sparse Matrix Decomposition via the Truncated Nuclear Norm and a Sparse Regularizer,"
% submitted to The Visual Computer, April 2018.
%
% Written by Zhichao Xue, version 1.0

clear all;clc;

%% load data
name='yale05';
X = load([name,'.mat']);
X = X.X;

%% process
sizem = size(X,1);
sizen = size(X,2);

iter = zeros(10,3);
err = zeros(10,3);
tim = zeros(10,3);

lower_R = 1; upper_R = 1;  
fprintf('now is running admm optimization method to recovery low-rank image and sparse image\n');

cd ADMM
tic;
[admmret]= admm_pic(X,lower_R,upper_R);
cd ..
%toc;
admm_num_iteration = max(admmret.iterations);
admm_time_cost = toc;

fprintf('\n TNNR-admm: time(%.1fs), iterations(%d)\n',admm_time_cost,admm_num_iteration);

