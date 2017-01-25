function [y] = functionPotentialWell(x,a,b,c)
    y = c+(x-a)* (x+a)*x*x*b;
end
