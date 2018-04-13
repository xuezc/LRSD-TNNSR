
function [ Z,E,iterations,LRerr,Sperr] = admmAXB(Z0, E0, A,B,X,M,rhomax,p)

% solve  minimize ||Z||_*-max(trace(A*Z*B'))+lamada*||E||+gamma*||W||

MAX_ITER = 200;

[m,n] = size(M);
lamda = 1/(max(m,n))^(1/2);
gamma = 0.09/(min(m,n))^(1/2);
rho = 0.25/mean(X(:));

Z=zeros(m,n);
E=zeros(m,n);
W=zeros(m,n);
Y=zeros(m,n);
P=zeros(m,n);


AB = A'*B;
eeppss = 1e-9;

for k = 1:MAX_ITER
    % Z-update
    lastZ = Z;
    tem1 = Y + AB;
    tem2 =  X - E + tem1/rho;
    tem3 = ifft2( W + P/rho );
    tem4 = 1/2 * ( tem2 + tem3 );
    [u,sigma,v] = svd(tem4);
    Z = u*max(sigma-1/rho/2,0)*v';
    N1 = norm(Z-lastZ,'fro');
    if(N1<eeppss)
        break;
    end
    
    %E-update
    lastE = E;
    tem5 =  X - Z + Y/rho;
    E = sign(tem5) .* max(abs(tem5)-lamda/rho,0);
    N2 = norm(E-lastE,'fro');
    if(N2<eeppss)
        break;
    end

    %W-update
    lastW = W;
    tem6 = fft2(Z);
    tem7 = tem6 - P/rho;
    W = sign(tem7) .* max(abs(tem7)-gamma/rho,0);
    N3 = norm(W-lastW,'fro');
    if(N3<eeppss)
        break;
    end   
 
    %Y-update
    Y = Y + rho * ( X - Z - E );
    
    %P-update
    P = P + rho * ( W - tem6 );
    
    % rho-update
    rho = min(p*rho,rhomax);
    
       
end

iterations = k;
LRerr = norm(Z-Z0,'fro')/norm(Z0,'fro');
Sperr = norm(E-E0,'fro')/norm(E0,'fro');

end
    
    
    
