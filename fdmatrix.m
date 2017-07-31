function [D] = fdmatrix(x,diff_ord,accu_ord)
%Function to output a finite difference matrix D based on the function
%fdcoefs
%x is the grid of nodes [x0 x1 x2...xN]
%diff_ord is the differentiation order
%accu_ord is the accuracy order

n = accu_ord+diff_ord-1;    %n+1 is the number of points in the stencil
m = diff_ord;               %m is the differentiation order

npoints = length(x);        %Number of grid points
nedge = floor((n+1)/2);     %Number of end points requiring consideration 

D = zeros(npoints);
if(mod(n,2)==1)
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(n,m,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge-1)) = fdcoefs(n,m,x((i-nedge):(i+nedge-1)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(n,m,x((end-n):end),x(i));
    end
else
    for i = 1:nedge
        D(i,1:(n+1)) = fdcoefs(n,m,x(1:(n+1)),x(i));
    end
    for i = (nedge+1):(npoints-nedge)
        D(i,(i-nedge):(i+nedge)) = fdcoefs(n,m,x((i-nedge):(i+nedge)),x(i));
    end
    for i = (npoints-nedge+1):(npoints)
        D(i,(end-n):end) = fdcoefs(n,m,x((end-n):end),x(i));
    end
end

D = sparse(D);
end