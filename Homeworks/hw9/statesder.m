function xd = statesder(t,z,u)
    x1 = z(1);
    x2 = z(2);
    
    x1d = -2*(x1+0.25) + (x2+0.5)*exp(25*x1/(x1+2));
    x2d = 0.5 - x2 - (x2+0.5)*exp(25*x1/(x1+2));
    
    xd = [x1d; x2d];
end