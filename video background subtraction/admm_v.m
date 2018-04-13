
function [Z,E] = admm_v(X,lower_R,upper_R)
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
% Z: the low-rank matrix
% E: the sparse matrix

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

Xfull = X;
[m, n] = size(Xfull);
h = min(m,n);
time_cost = zeros(h,1);
iterations_cost = zeros(h,1);

number_of_out_iter = 1;

fprintf('now is running under rank(r)=                   ');


for R = lower_R:upper_R
    tic;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%2d',R);
    X = Xfull(:,:);
    M = X;
    X_rec = zeros(size(Xfull));
    Z_rec = zeros(size(Xfull));
    E_rec = zeros(size(Xfull));
        
    for out_iter = 1:number_of_out_iter
        %out_iter
        [u, sigma, v] = svd(X,'econ');
         A = u(:,1:R)'; B = v(:,1:R)';
            
        [Z,E,W,iter_count] = admmAXB(A,B,X,M,0.01,1.1);
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

Z = Z ./ 255;

E = E ./ 255;

end