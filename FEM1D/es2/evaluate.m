function u0 = evaluate(u_h,x_h)
    
    % define the point of the domain where to evaluate the approx. solution
    x0 = 1;
    idx = find(x_h==x0);
    u0 = u_h(idx);
    
end