function [gradF] = gradientPerturbation(x,H,gamma)
    gradF = H * x - gamma.* sin(x);
end
