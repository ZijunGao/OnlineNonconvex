% perturbed quadratic objective function

% objective function (functionPerturbation) : f(x) = x' * H * x / 2 +  gamma * (cos(x(1)) ... + cos(x(d))), where H
% is a d by d matrix, Hessian f has negative minimal eigenvalue, and x is a d
% dimensional vector.

% SGD iterates : (1) sample independent  Gaussian noise/ Bernoulli noise/ Uniform noise/ Uniform-spherical noise with total variance sigma^2; 
% (2) x  = x - grad f(x) + eta * noise. 
% X_SGD = X_SGD - eta * grad f(X_SGD) + sigma * eta *
% randn(d, 1)/ 2 * sqrt(3) * (rand(d)-0.5) / sign(randn(d)) /  sqrt(d) *
% randn(d, 1)/||rand(d, 1)||_2.

% stopping rule: function value decreases below a certain level.

% plot of log(iteration number) - log(1/eta).

%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H = diag([10,10,10,10,-10]); % base of objective function

s = 4; % number of gammas
gamma = 10./2.^(0:(s-1))-10;% magnitude of perturbation
n = 10; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -1 + functionPerturbation(x0, H, gamma); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(n)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + eta(n) * abs(H(d,d)));
 N = floor (min(1.5*1e5, 10 * NL));

m = 200; % simulation times

%% iteration

 type = 'Gaussian';
% type = 'Uniform';
% type = 'Bernoulli';
% type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0; % the same initialization

for t = 1 : s
    functionLevel = -1 + functionPerturbation(x0, H, gamma(t)); % objective value of function
    for l = 1 : n
        for k = 1 : m
            functionValue = functionPerturbation(X_SGD(:, 1), H, gamma(t));
            i = 1;
            while functionValue > functionLevel
               i = i + 1;
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPerturbation(X_SGD(:, i - 1), H, gamma(t)) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % Uniform noise  
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % Spherical Gaussian noise
               functionValue = functionPerturbation(X_SGD(:, i), H, gamma(t));
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
%            disp(k);
        end
%       disp(recordIterationNumber(t, l,:));
%       disp(mean(recordIterationNumber(t, l,:))); 
%       disp(std(recordIterationNumber(t, l,:))); 
%       hist(recordIterationNumber(t, l,:));
        disp(l);
    end
    disp(t);
end

%% history
 save perturbationGauDat;
% save perturbationUniDat;
% save perturbationBerDat;
% save perturbationSphDat;



















%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H = diag([10,10,10,10,-10]); % base of objective function

s = 4; % number of gammas
gamma = 10./2.^(0:(s-1))-10;% magnitude of perturbation
n = 10; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -1 + functionPerturbation(x0, H, gamma); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(n)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + eta(n) * abs(H(d,d)));
 N = floor (min(1.5*1e5, 10 * NL));

m = 200; % simulation times

%% iteration

% type = 'Gaussian';
 type = 'Uniform';
% type = 'Bernoulli';
% type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0; % the same initialization

for t = 1 : s
    functionLevel = -1 + functionPerturbation(x0, H, gamma(t)); % objective value of function
    for l = 1 : n
        for k = 1 : m
            functionValue = functionPerturbation(X_SGD(:, 1), H, gamma(t));
            i = 1;
            while functionValue > functionLevel
               i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPerturbation(X_SGD(:, i - 1), H, gamma(t)) + eta(l) * sigma * randn(d, 1); % Gaussian noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % Uniform noise  
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % Spherical Gaussian noise
               functionValue = functionPerturbation(X_SGD(:, i), H, gamma(t));
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
%            disp(k);
        end
%       disp(recordIterationNumber(t, l,:));
%       disp(mean(recordIterationNumber(t, l,:))); 
%       disp(std(recordIterationNumber(t, l,:))); 
%       hist(recordIterationNumber(t, l,:));
        disp(l);
    end
    disp(t);
end

%% history
% save perturbationGauDat;
 save perturbationUniDat;
% save perturbationBerDat;
% save perturbationSphDat;




















%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H = diag([10,10,10,10,-10]); % base of objective function

s = 4; % number of gammas
gamma = 10./2.^(0:(s-1))-10;% magnitude of perturbation
n = 10; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -1 + functionPerturbation(x0, H, gamma); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(n)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + eta(n) * abs(H(d,d)));
 N = floor (min(1.5*1e5, 10 * NL));

m = 200; % simulation times

%% iteration

% type = 'Gaussian';
% type = 'Uniform';
 type = 'Bernoulli';
% type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0; % the same initialization

for t = 1 : s
    functionLevel = -1 + functionPerturbation(x0, H, gamma(t)); % objective value of function
    for l = 1 : n
        for k = 1 : m
            functionValue = functionPerturbation(X_SGD(:, 1), H, gamma(t));
            i = 1;
            while functionValue > functionLevel
               i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPerturbation(X_SGD(:, i - 1), H, gamma(t)) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % Uniform noise  
               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % Spherical Gaussian noise
               functionValue = functionPerturbation(X_SGD(:, i), H, gamma(t));
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
%            disp(k);
        end
%       disp(recordIterationNumber(t, l,:));
%       disp(mean(recordIterationNumber(t, l,:))); 
%       disp(std(recordIterationNumber(t, l,:))); 
%       hist(recordIterationNumber(t, l,:));
        disp(l);
    end
    disp(t);
end

%% history
% save perturbationGauDat;
% save perturbationUniDat;
 save perturbationBerDat;
% save perturbationSphDat;





















%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H = diag([10,10,10,10,-10]); % base of objective function

s = 4; % number of gammas
gamma = 10./2.^(0:(s-1))-10;% magnitude of perturbation
n = 10; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -1 + functionPerturbation(x0, H, gamma); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(n)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + eta(n) * abs(H(d,d)));
 N = floor (min(1.5*1e5, 10 * NL));

m = 200; % simulation times

%% iteration

% type = 'Gaussian';
% type = 'Uniform';
% type = 'Bernoulli';
 type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0; % the same initialization

for t = 1 : s
    functionLevel = -1 + functionPerturbation(x0, H, gamma(t)); % objective value of function
    for l = 1 : n
        for k = 1 : m
            functionValue = functionPerturbation(X_SGD(:, 1), H, gamma(t));
            i = 1;
            while functionValue > functionLevel
               i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPerturbation(X_SGD(:, i - 1), H, gamma(t)) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % Uniform noise  
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sign(randn(d, 1)); % Bernoulli noise
               noiseSph = randn(d, 1);
               X_SGD(:, i) = X_SGD(:, i - 1) - eta * gradientPerturbation(X_SGD(:, i - 1), H, gamma) + eta * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % Spherical Gaussian noise
               functionValue = functionPerturbation(X_SGD(:, i), H, gamma(t));
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
%            disp(k);
        end
%       disp(recordIterationNumber(t, l,:));
%       disp(mean(recordIterationNumber(t, l,:))); 
%       disp(std(recordIterationNumber(t, l,:))); 
%       hist(recordIterationNumber(t, l,:));
        disp(l);
    end
    disp(t);
end

%% history
% save perturbationGauDat;
% save perturbationUniDat;
% save perturbationBerDat;
 save perturbationSphDat;

%% bugs and tips

