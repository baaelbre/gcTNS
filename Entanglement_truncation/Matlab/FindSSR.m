function [rhor,lambda,flag]=FindSSR(Q,R,D,maxiter,tol,rho0)

if nargin<4 || isempty(maxiter)
    maxiter=200;
end
if nargin<5 || isempty(tol)
    tol=eps;
end

% Right action of transfer operator
    function w=ApplyTR(v)
        rho=Vec2SMat(D,v);
        rho=Q*rho+rho*Q'+R*rho*R';
        rho=(rho+rho')/2;
        w=SMat2Vec(D,rho);
    end

% Compute eigenvalue and eigenvector
    opts.disp=0;
    opts.tol=max(tol/100,eps);
    opts.maxit=maxiter;
    opts.isreal=true;
    opts.issym=false;
    if nargin>=6 && ~isempty(rho0)
        opts.v0=SMat2Vec(D,rho0);
    end
    
    try
        [vr,lambda,flag]=eigs(@ApplyTR,D*(D+1)/2,[],1,'LR',opts);
    catch err
        opts.tol=max(tol,eps);
        opts.maxit=10*maxiter;
        warning('Trying less accurate determination of right eigenvector');
        [vr,lambda,flag]=eigs(@ApplyTR,D*(D+1)/2,[],1,'LR',opts);
    end
    rhor=Vec2SMat(D,vr);
    rhor=rhor/trace(rhor);
end

% Functions to convert symmetric matrix into vector and vice versa
function M=Vec2SMat(D,V)
if length(V)==D*(D+1)/2
    ind=1;
    for i=1:D
        M(i,i)=V(ind);
        ind=ind+1;
        for j=i+1:D
            M(i,j)=V(ind)/sqrt(2);
            M(j,i)=V(ind)/sqrt(2);
            ind=ind+1;
        end
    end
end
end

function V=SMat2Vec(D,M)
V=zeros(D*(D+1)/2,1);
ind=1;
for i=1:D
    V(ind)=M(i,i);
    ind=ind+1;
    for j=i+1:D
        V(ind)=(M(i,j)+M(j,i))/sqrt(2);
        ind=ind+1;
    end
end
end