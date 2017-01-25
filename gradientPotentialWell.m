function [y] = gradientPotentialWell(x,a,b,c)
    y = b*((x+a)*x*x+(x-a)*x*x+2*(x+a)*(x-a)*x);
end
