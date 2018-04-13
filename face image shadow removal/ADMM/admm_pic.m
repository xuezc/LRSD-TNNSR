
function [ret] = admm_pic(X,lower_R,upper_R)
% This function implements the LRSD-TNNSR algorithm 
% @ April 2018, Nanjing Tech University
% 
% Details about the algorithm can be found in the paper:
% Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali,
% "Low-rank and Sparse Matrix Decomposition via the Truncated Nuclear Norm and a Sparse Regularizer,"
% submitted to The Visual Computer, April 2018.
%
% Inputs:
% X: the data matrix
% lower_R: lower rank
% upper_R: higher rank

% Outputs:
% ret.time: the running time
% ret.iterations: iteration number

% Written by Zhichao Xue,version 1.0                                    
%
% If you have any questions or comments regarding this package, or if you want to 
% report any bugs or unexpected error messages, please send an e-mail to
% xuezhichao@njtech.edu.cn
%     
% Copyright 2018 Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali
% 
% This software is a free software distributed under the terms of the GNU 
% Public License version 3 (http://www.gnu.org/licenses/gpl.txt). You can 
% redistribute it and/or modify it under the terms of this licence, for 
% personal and non-commercial use and research purpose.

[m, n] = size(X);
h = min(m,n);
time_cost = zeros(h,1);
iterations_cost = zeros(h,1);

number_of_out_iter = 1;

fprintf('now is running under rank(r)=                   ');

X1 = X(:,30);
X2 = reshape(X1,192,168);

figure(1);
subplot(1,1,1);
imshow(X2,[]);
xlabel('original image');

for i=1:64
    X1i = X(:,i);
    X2i = reshape(X1i,192,168);
    imwrite(im2uint8(mat2gray(X2i)),['input\','X',num2str(i),'.jpg']);
end

for R = lower_R:upper_R
    tic;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%2d',R);
    X = X(:,:);
    
    X_rec = zeros(size(X));
    Z_rec = zeros(size(X));
    E_rec = zeros(size(X));
        
    for out_iter = 1:number_of_out_iter
        %out_iter
        [u, sigma, v] = svd(X,'econ');
         A = u(:,1:R)'; B = v(:,1:R)';
            
        [Z,E,W,iter_count] = admmAXB(A,B,X,10000,1.01);
         lastX = X;
         X = Z + E;

         Z_rec(:,:,out_iter) = Z;
         E_rec(:,:,out_iter) = E;
         X_rec(:,:,out_iter) = X;
         
         iterations_cost(R) = iterations_cost(R) + iter_count;

         if(out_iter>=2 && norm(X_rec(:,:,out_iter)-X_rec(:,:,out_iter-1),'fro')<1e-7)
             X = X_rec(:,:,out_iter);
             break;
         end
         X = X_rec(:,:,out_iter);
         Z = Z_rec(:,:,out_iter);
         E = E_rec(:,:,out_iter);
    end
    
    time_cost(R) = toc;

end

Z1 = Z(:,30);
Z2 = reshape(Z1,192,168);
E1 = E(:,30);
E2 = reshape(E1,192,168);
imwrite(im2uint8(mat2gray(Z2)),['output\','Z.jpg']);
imwrite(im2uint8(mat2gray(E2)),['output\','E.jpg']);

figure(2);
subplot(1,1,1); % show the recovery low-rank image
imshow(Z2,[]);
xlabel('recovered low-rank image');

figure(3);
subplot(1,1,1); % show the recovery sparse image
imshow(E2,[]);
xlabel('recovered sparse image');

for i=1:64
    Z1i = Z(:,i);
    Z2i = reshape(Z1i,192,168);
    imwrite(im2uint8(mat2gray(Z2i)),['output\','Z',num2str(i),'.jpg']);
end

for i=1:64
    E1i = E(:,i);
    E2i = reshape(E1i,192,168);
    imwrite(im2uint8(mat2gray(E2i)),['output\','E',num2str(i),'.jpg']);
end

ret.time = time_cost;
ret.iterations = iterations_cost;
end