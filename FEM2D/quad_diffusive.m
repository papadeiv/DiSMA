function sum = quad_diffusive(grad_j,grad_k,B)
    
    global coefficient_functions;
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        sum = sum + omega(q)*coefficient_functions{1}(xhat(q),yhat(q))...
                            *transpose([grad_j{1}(xhat(q),yhat(q));grad_j{2}(xhat(q),yhat(q))])...
                            *B*[grad_k{1}(xhat(q),yhat(q));grad_k{2}(xhat(q),yhat(q))];
    end

end