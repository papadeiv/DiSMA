function [d_x,d_y] = delta(nodes, triangle)

    % assign the e-th triangle's local vertices coordinates to v1,v2,v3
    v1 = nodes(triangle(1),3:4);
    v2 = nodes(triangle(2),3:4);
    v3 = nodes(triangle(3),3:4);

    x = [v1(1) v2(1) v3(1)];
    y = [v1(2) v2(2) v3(2)];
    
    d_x(1) = x(3)-x(2);
    d_x(2) = x(1)-x(3);
    d_x(3) = x(2)-x(1);
    d_y(1) = y(2)-y(3);
    d_y(2) = y(3)-y(1);
    d_y(3) = y(1)-y(2);
    
end

