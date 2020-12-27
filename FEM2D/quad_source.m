function sum = quad_source(phi)

    global source_function;
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        sum = sum + omega(q)*source_function(xhat(q),yhat(q))*phi(xhat(q),yhat(q));
    end
    

end