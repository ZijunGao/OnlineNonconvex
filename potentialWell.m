% potential well problem

% objective function (functionPotentialWell) : f(x) = b * (x^2 - a^2) * x^2
% + c.

% SGD iterates : (1) sample independent  Gaussian noise/ Bernoulli noise/ Uniform noise/ Uniform-spherical noise with total variance sigma^2; 
% (2) x  = x - eta * H * x + eta * noise. 
% X_SGD = X_SGD - eta * H * X_SGD + sigma * eta *
% randn(d, 1)/ 2 * sqrt(3) * (rand(d)-0.5) / sign(randn(d)) /  sqrt(d) *
% randn(d, 1)/||rand(d, 1)||_2.

% stopping rule: function value decreases below a certain level.

%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 1; % dimension
s = 4; % number of minimal eigenvalues, or b
a = 1;
b = 2./2.^(0:(s-1));% magnitude of minimal eigenvalue
c = 0;
n = 8; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -0.1 * a^4 * b + functionPotentialWell(x0,a,b,c); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + 2*a^2*b(1)*eta(1));
 N = floor (min(1.5*1e5, 1.2 * NL));

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
    for l = 1 : n
        functionLevel = -0.25 * a^4 * b(t) + functionPotentialWell(x0,a,b(t),c) + eta(l); % objective value of function 
        for k = 1 : m
            functionValue = functionPotentialWell(X_SGD(:, 1), a,b(t),c);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionPotentialWell(X_SGD(:, i), a,b(t),c);
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
            disp(k);
        end
        disp(l);
%       disp(recordIterationNumber(l,:));
%       disp(mean(recordIterationNumber(l,:))); 
%       disp(std(recordIterationNumber(l,:))); 
%       hist(recordIterationNumber(l,:));
    end
    disp(t);
end
disp(type);

%% history
 save potentialWellGauDat;
% save potentialWellUniDat;
% save potentialWellBerDat;
% save potentialWellSphDat;


















%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 1; % dimension
s = 4; % number of minimal eigenvalues, or b
a = 1;
b = 2./2.^(0:(s-1));% magnitude of minimal eigenvalue
c = 0;
n = 8; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -0.1 * a^4 * b + functionPotentialWell(x0,a,b,c); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + 2*a^2*b(1)*eta(1));
 N = floor (min(1.5*1e5, 1.2 * NL));

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
    for l = 1 : n
        functionLevel = -0.25 * a^4 * b(t) + functionPotentialWell(x0,a,b(t),c) + eta(l); % objective value of function 
        for k = 1 : m
            functionValue = functionPotentialWell(X_SGD(:, 1), a,b(t),c);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * randn(d, 1); % Gaussian noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionPotentialWell(X_SGD(:, i), a,b(t),c);
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
            disp(k);
        end
        disp(l);
%       disp(recordIterationNumber(l,:));
%       disp(mean(recordIterationNumber(l,:))); 
%       disp(std(recordIterationNumber(l,:))); 
%       hist(recordIterationNumber(l,:));
    end
    disp(t);
end
disp(type);

%% history
% save potentialWellGauDat;
 save potentialWellUniDat;
% save potentialWellBerDat;
% save potentialWellSphDat;


















%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 1; % dimension
s = 4; % number of minimal eigenvalues, or b
a = 1;
b = 2./2.^(0:(s-1));% magnitude of minimal eigenvalue
c = 0;
n = 8; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -0.1 * a^4 * b + functionPotentialWell(x0,a,b,c); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + 2*a^2*b(1)*eta(1));
 N = floor (min(1.5*1e5, 1.2 * NL));

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
    for l = 1 : n
        functionLevel = -0.25 * a^4 * b(t) + functionPotentialWell(x0,a,b(t),c) + eta(l); % objective value of function 
        for k = 1 : m
            functionValue = functionPotentialWell(X_SGD(:, 1), a,b(t),c);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionPotentialWell(X_SGD(:, i), a,b(t),c);
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
            disp(k);
        end
        disp(l);
%       disp(recordIterationNumber(l,:));
%       disp(mean(recordIterationNumber(l,:))); 
%       disp(std(recordIterationNumber(l,:))); 
%       hist(recordIterationNumber(l,:));
    end
    disp(t);
end
disp(type);

%% history
% save potentialWellGauDat;
% save potentialWellUniDat;
 save potentialWellBerDat;
% save potentialWellSphDat;


























%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 1; % dimension
s = 4; % number of minimal eigenvalues, or b
a = 1;
b = 2./2.^(0:(s-1));% magnitude of minimal eigenvalue
c = 0;
n = 8; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
% functionLevel = -0.1 * a^4 * b + functionPotentialWell(x0,a,b,c); % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(1) / sigma^2)) / log (1 + 2*a^2*b(1)*eta(1));
 N = floor (min(1.5*1e5, 1.2 * NL));

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
    for l = 1 : n
        functionLevel = -0.25 * a^4 * b(t) + functionPotentialWell(x0,a,b(t),c) + eta(l); % objective value of function 
        for k = 1 : m
            functionValue = functionPotentialWell(X_SGD(:, 1), a,b(t),c);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
               noiseSph = randn(d, 1);
               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * gradientPotentialWell(X_SGD(:, i - 1), a,b(t),c) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionPotentialWell(X_SGD(:, i), a,b(t),c);
%               disp(functionValue)
            end
            recordIterationNumber(t, l, k) = i; 
            disp(k);
        end
        disp(l);
%       disp(recordIterationNumber(l,:));
%       disp(mean(recordIterationNumber(l,:))); 
%       disp(std(recordIterationNumber(l,:))); 
%       hist(recordIterationNumber(l,:));
    end
    disp(t);
end
disp(type);

%% history
% save potentialWellGauDat;
% save potentialWellUniDat;
% save potentialWellBerDat;
 save potentialWellSphDat;

%% bugs and tips



