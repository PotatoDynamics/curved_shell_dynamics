function [dY] = EOM(tau,x,f_N,ome_N,n,dUs,Q)

    xi = x(1:n/2);
    dxi = x(n/2+1:end);
    
    y = dxi;
    dy = Q(dxi,f_N,ome_N,tau)-dUs(xi,tau);

    dY = [y; dy];