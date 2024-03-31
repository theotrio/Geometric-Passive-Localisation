function [x,xlist,iter] = Newton_Raphson_Mod(f,df,initp,tol,maxit,omega,xd,yd,beta1,xb,yb,x0,y0,phi,thetam)
    if nargin<3
        error('Not enough data\n');
    end
    if nargin<4
        tol = eps;
    end
    if nargin<5
        maxit = 50;
    end
    iter = 0;
    x=initp;
    xlist = x;
    xdiff = inf;
    while xdiff>=tol
        iter=iter+1;
        xold = x;
        x = x - feval(f,omega,xd,yd,beta1,xb,yb,x0,y0,phi,thetam,x)/feval(df,omega,xd,yd,beta1,xb,yb,x0,y0,phi,thetam,x);
        xdiff = abs(x-xold)/abs(x);
        xlist=[xlist;x];
        if iter>=maxit
            x=-inf;
            %error('Not converged after max iterations.');
        end
    end
end