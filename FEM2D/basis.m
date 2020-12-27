function [phi, grad_phi] = basis()

    phi = {
        @(x,y) (1-x-y)*(2*(1-x-y)-1);
        @(x,y) x*(2*x-1);
        @(x,y) y*(2*y-1);
        @(x,y) 4*(1-x-y)*x;
        @(x,y) 4*x*y;
        @(x,y) 4*y*(1-x-y);
        };
    
    grad_phi = {
        @(x,y) -3+4*x+4*y, @(x,y) -3+4*y+4*x;
        @(x,y) 4*x-1, @(x,y) 0;
        @(x,y) 0, @(x,y) 4*y-1;
        @(x,y) 4*(1-2*x-y), @(x,y) -4*x;
        @(x,y) 4*y, @(x,y) 4*x;
        @(x,y) -4*y, @(x,y) 4*(1-2*y-x);
        };

end