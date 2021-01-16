function sum = quad_source(j, timestep)
    
    global source_function;
    [phihat, grad_phihat] = basis();
    phi = phihat{j};
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        sum = sum + omega(q)*source_function(xhat(q),yhat(q),timestep)*phi(xhat(q),yhat(q));
    end
    
end