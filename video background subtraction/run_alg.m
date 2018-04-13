lower_R = 1; upper_R = 1;  
fprintf('now is running admm optimization method to recovery low-rank image and sparse image\n');

[Z,E] = admm_v(M,lower_R,upper_R);

L = Z;
S = E; 