function sum = quad_source(j)
    
    global t;
    global source_function;
    [phihat, grad_phihat] = basis();
    phi = phihat{j};
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        sum = sum + omega(q)*source_function(xhat(q),yhat(q),t)*phi(xhat(q),yhat(q));
    end
    
end