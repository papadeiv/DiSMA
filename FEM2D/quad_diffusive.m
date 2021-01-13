function sum = quad_diffusive(j,k,diffB)
    
    global t;
    global coefficient_functions;
    % define basis functions on the reference triangle
    [phihat, grad_phihat] = basis();
    % import nodes and weigths on the reference traingle
    [xhat,yhat,omega] = quadrature_weights;
    Nq = length(xhat);
    sum = 0;
    for q=1:Nq
        sum = sum + omega(q)*coefficient_functions{1}(xhat(q),yhat(q),t)*...
                            [grad_phihat{j,1}(xhat(q),yhat(q)),grad_phihat{j,2}(xhat(q),yhat(q))]*...
                            diffB*...
                            [grad_phihat{k,1}(xhat(q),yhat(q)),grad_phihat{k,2}(xhat(q),yhat(q))]';
    end

end