% The demo of synthetic data

% References:
% Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali,
% "Low-rank and Sparse Matrix Decomposition via the Truncated Nuclear Norm and a Sparse Regularizer,"
% submitted to The Visual Computer, April 2018.
%
% Written by Zhichao Xue, version 1.0

clear all;clc;

%% load synthetic data
name='100';
A = load([name,'.mat']);
X=A.X;
Z=A.Z;
E=A.E;

%% show image
figure(1);
subplot(1,1,1);
imshow(X,[]);
xlabel('original image');
imwrite(im2uint8(mat2gray(X)),['save\','X0.jpg']);    
figure(2);
subplot(1,1,1);
imshow(Z,[]);
xlabel('Low-Rank image');
imwrite(im2uint8(mat2gray(Z)),['save\','Z0.jpg']);
figure(3);
subplot(1,1,1);
imshow(E,[]);
xlabel('Sparse image');
imwrite(im2uint8(mat2gray(E)),['save\','E0.jpg']);

%% process
sizem = size(X,1);
sizen = size(X,2);

iter = zeros(10,3);
err = zeros(10,3);
tim = zeros(10,3);

lower_R = 10; upper_R = 10;  
fprintf('now is running admm optimization method to recovery low-rank image and sparse image\n');

Z0 = Z;
E0 = E;

cd ADMM
tic;
[admmret]= admm_pic(Z0, E0, X,lower_R,upper_R);
cd ..
%toc;
admm_num_iteration = max(admmret.iterations);
admm_time_cost = toc;

fprintf('\n TNNR-admm: time(%.4fs),Totalerr(%.15f),LRerr(%.15f),Sperr(%.15f), iterations(%d)\n',admm_time_cost,admmret.Totalerr,admmret.LRerr,admmret.Sperr,admm_num_iteration);

