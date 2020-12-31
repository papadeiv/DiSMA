function sum = quad_source(j,dx,dy,v1)

    global source_function;
    [phihat, grad_phihat] = basis();
    phi = phihat{j};
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        x_mapped = v1(1) + dx(3)*xhat(q) - dx(2)*yhat(q);
        y_mapped = v1(2) + dy(2)*yhat(q) - dy(3)*xhat(q);
        sum = sum + omega(q)*source_function(x_mapped,y_mapped)*phi(xhat(q),yhat(q));
    end
    
end