
function [ Z,E,W,iterations] = admmAXB( A,B,X,rhomax,p)

% solve  minimize ||Z||_*-max(trace(A*Z*B'))+lamada*||E||_1+gamma*||W||_1

MAX_ITER = 40;

[m,n] = size(X);
lamda = 0.005/(min(m,n))^(1/2);
gamma = 0.0009/(min(m,n))^(1/2);
rho = 0.005;             

Z=zeros(m,n);
E=zeros(m,n);
W=zeros(m,n);
Y=zeros(m,n);
P=zeros(m,n);


AB = A'*B;
eeppss = 1e-7;

for k = 1:MAX_ITER
    % Z-update
    lastZ = Z;
    tem1 = Y + AB;
    tem2 =  X - E + tem1/rho;
    tem3 = mirt_idctn( W + P/rho );
    tem4 = 1/2 * ( tem2 + tem3 );
    [u,sigma,v] = svd(tem4,'econ');
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
    tem6 = mirt_dctn(Z);
    tem7 = tem6 - P/rho;
    W = sign(tem7) .* max(abs(tem7)-gamma/rho,0);
    N2 = norm(W-lastW,'fro');
    if(N2<eeppss)
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

end
    
    
    
