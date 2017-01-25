function [y] = functionPerturbation(x,H,gamma)
    y = x' * H * x /2 + gamma * sum(cos(x));
end
