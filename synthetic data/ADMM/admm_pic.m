
function [ret, E] = admm_pic(Z0, E0, X, lower_R, upper_R)
% This function implements the LRSD-TNNSR algorithm 
% @ April 2018, Nanjing Tech University
% 
% Details about the algorithm can be found in the paper:
% Z. Xue, J. Dong, Y. Zhao, C. Liu, and R. Chellali,
% "Low-rank and Sparse Matrix Decomposition via the Truncated Nuclear Norm and a Sparse Regularizer,"
% submitted to The Visual Computer, April 2018.
%
% Inputs:
% Z0: original low-rank matrix
% E0: original sparse matrix
% X: original complete matrix
% lower_R: lower rank
% upper_R: higher rank

% Outputs:
% ret.time: the running time
% ret.iterations: iteration number
% ret.LRerr: the error of the low-rank component
% ret.Sperr: the error of the sparse component
% ret.Totalerr: the error of the complete matrix

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

Zrecover = zeros(m,n,h);
Erecover = zeros(m,n,h);
number_of_out_iter = 1;

fprintf('now is running under rank(r)=                   ');

for R = lower_R:upper_R
    tic;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%2d',R);
   
    M = X;
        
    for out_iter = 1:number_of_out_iter
        %out_iter
        [u, sigma, v] = svd(X);
         A = u(:,1:R)'; B = v(:,1:R)';
            
        [Z,E,iter_count,LRerr,Sperr] = admmAXB(Z0, E0,A,B,X,M,0.01,1.1);
         lastX = X;
         X = Z + E;       
         Z_rec(:,:,out_iter) = Z;
         E_rec(:,:,out_iter) = E;
         X_rec(:,:,out_iter) = X;
        
         iterations_cost(R) = iterations_cost(R) + iter_count;        
         
         X = X_rec(:,:,out_iter);
         Z = Z_rec(:,:,out_iter);
         E = E_rec(:,:,out_iter);
         Totalerr = norm(X-lastX,'fro')/norm(lastX,'fro');
         if(Totalerr<1e-9)
             break
         end
         
    end   

    
    time_cost(R) = toc;    
end
figure(4);
subplot(1,1,1); % show the recovery low-rank image
Zrecover=max(Zrecover(:,:,R),0);
Zrecover = min(Zrecover,255);
imshow(Z./255,[]);
xlabel('recovered low-rank image');
imwrite(im2uint8(mat2gray(Z)),['result\','Z.jpg']);

figure(5);
subplot(1,1,1); % show the recovery sparse image
Erecover=max(Erecover(:,:,R),0);
Erecover = min(Erecover,255);
imshow(E./255,[]);
xlabel('recovered sparse image');
imwrite(im2uint8(mat2gray(E)),['result\','E.jpg']);

ret.time = time_cost;
ret.iterations = iterations_cost;
ret.LRerr = LRerr;
ret.Sperr = Sperr;
ret.Totalerr = Totalerr;
end