% quadratic objective function

% objective function (functionQuadratic) : f(x) = x' * H * x / 2, where H
% is a d by d matrix and has negative minimal eigenvalue, and x is a d
% dimensional vector.

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

d = 5; % dimension
H0 = diag([1,1,1,1,-1]); % base of objective function
s = 4; % number of minimal eigenvalues
gamma = 10./2.^(0:(s-1));% magnitude of minimal eigenvalue
n = 10; % number of eta
eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
functionLevel = -1; % objective value of function 

 alpha = 0.1; % tolerance 
 NL = ((0.5 + alpha) * log(1/eta(n)) + 0.5 * log(2 * abs(functionLevel) / sigma^2)) / log (1 + eta(n) * abs(H0(d,d)));
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
    H  =  H0.* gamma(t);
    for l = 1 : n
        for k = 1 : m
            functionValue = functionQuadratic(X_SGD(:, 1), H);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionQuadratic(X_SGD(:, i), H);
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
 save quadraticGauDat;
% save quadraticUniDat;
% save quadraticBerDat;
% save quadraticSphDat;

%% bugs and tips
% 1. any bug is possible - be brave and patient to debug
















% quadratic objective function

% objective function (functionQuadratic) : f(x) = x' * H * x / 2, where H
% is a d by d matrix and has negative minimal eigenvalue, and x is a d
% dimensional vector.

% SGD iterates : (1) sample independent  Gaussian noise/ Bernoulli noise/ Uniform noise/ Uniform-spherical noise with total variance sigma^2; 
% (2) x  = x - eta * H * x + eta * noise. 
% X_SGD = X_SGD - eta * H * X_SGD + sigma * eta *
% randn(d, 1)/ 2 * sqrt(3) * (rand(d)-0.5) / sign(randn(d)) /  sqrt(d) *
% randn(d, 1)/||rand(d, 1)||_2.

% stopping rule: function value decreases below a certain level.

% plot of log(iteration number) - log(base/eta), where base is the largest
% eta used.

%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H0 = diag([1,1,1,1,-1]); % base of objective function
s = 4; % number of minimal eigenvalues
gamma = 10./2.^(0:(s-1));% magnitude of minimal eigenvalue
n = 10; % number of eta

eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
functionLevel = -1; % objective value of function 

alpha = 0.1; % tolerance 
NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(functionLevel) / sigma^2)) / log (1 + eta(1) * abs(H0(d,d)));
N = floor (min(1.5*1e4, 1.2 * NL));

m = 200; % simulation times

%% iteration

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0;

% type = 'Gaussian';
 type = 'Uniform';
% type = 'Bernoulli';
% type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);
for t = 1 : s
    H  =  H0.* gamma(t);
    for l = 1 : n
        for k = 1 : m
            functionValue = functionQuadratic(X_SGD(:, 1), H);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * randn(d, 1); % Gaussian noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionQuadratic(X_SGD(:, i), H);
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
% save quadraticGauDat;
 save quadraticUniDat;
% save quadraticBerDat;
% save quadraticSphDat;

%% bugs and tips
% 1. any bug is possible - be brave and patient to debug


















% quadratic objective function

% objective function (functionQuadratic) : f(x) = x' * H * x / 2, where H
% is a d by d matrix and has negative minimal eigenvalue, and x is a d
% dimensional vector.

% SGD iterates : (1) sample independent  Gaussian noise/ Bernoulli noise/ Uniform noise/ Uniform-spherical noise with total variance sigma^2; 
% (2) x  = x - eta * H * x + eta * noise. 
% X_SGD = X_SGD - eta * H * X_SGD + sigma * eta *
% randn(d, 1)/ 2 * sqrt(3) * (rand(d)-0.5) / sign(randn(d)) /  sqrt(d) *
% randn(d, 1)/||rand(d, 1)||_2.

% stopping rule: function value decreases below a certain level.

% plot of log(iteration number) - log(base/eta), where base is the largest
% eta used.

%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H0 = diag([1,1,1,1,-1]); % base of objective function
s = 4; % number of minimal eigenvalues
gamma = 10./2.^(0:(s-1));% magnitude of minimal eigenvalue
n = 10; % number of eta

eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
functionLevel = -1; % objective value of function 

alpha = 0.1; % tolerance 
NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(functionLevel) / sigma^2)) / log (1 + eta(1) * abs(H0(d,d)));
N = floor (min(1.5*1e4, 1.2 * NL));

m = 200; % simulation times

%% iteration

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0;

% type = 'Gaussian';
% type = 'Uniform';
 type = 'Bernoulli';
% type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);
for t = 1 : s
    H  =  H0.* gamma(t);
    for l = 1 : n
        for k = 1 : m
            functionValue = functionQuadratic(X_SGD(:, 1), H);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionQuadratic(X_SGD(:, i), H);
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
% save quadraticGauDat;
% save quadraticUniDat;
 save quadraticBerDat;
% save quadraticSphDat;

%% bugs and tips
% 1. any bug is possible - be brave and patient to debug




















% quadratic objective function

% objective function (functionQuadratic) : f(x) = x' * H * x / 2, where H
% is a d by d matrix and has negative minimal eigenvalue, and x is a d
% dimensional vector.

% SGD iterates : (1) sample independent  Gaussian noise/ Bernoulli noise/ Uniform noise/ Uniform-spherical noise with total variance sigma^2; 
% (2) x  = x - eta * H * x + eta * noise. 
% X_SGD = X_SGD - eta * H * X_SGD + sigma * eta *
% randn(d, 1)/ 2 * sqrt(3) * (rand(d)-0.5) / sign(randn(d)) /  sqrt(d) *
% randn(d, 1)/||rand(d, 1)||_2.

% stopping rule: function value decreases below a certain level.

% plot of log(iteration number) - log(base/eta), where base is the largest
% eta used.

%%
clear all; close all; clc;

%% parameter
seedn = 1000 * randn(1);
randn ('seed', seedn);
rand('seed', seedn);

d = 5; % dimension
H0 = diag([1,1,1,1,-1]); % base of objective function
s = 4; % number of minimal eigenvalues
gamma = 10./2.^(0:(s-1));% magnitude of minimal eigenvalue
n = 10; % number of eta

eta = 0.01./ 2.^(1:n); % learning rate
sigma = 1/sqrt(d); % variance of noise

x0 = zeros(d, 1); % initialization
functionLevel = -1; % objective value of function 

alpha = 0.1; % tolerance 
NL = ((0.5 + alpha) * log(1/eta(1)) + 0.5 * log(2 * abs(functionLevel) / sigma^2)) / log (1 + eta(1) * abs(H0(d,d)));
N = floor (min(1.5*1e4, 1.2 * NL));

m = 200; % simulation times

%% iteration

X_SGD = zeros(d , N);
X_SGD(:, 1) = x0;

% type = 'Gaussian';
% type = 'Uniform';
% type = 'Bernoulli';
 type = 'Uniform-spherical';

recordIterationNumber = zeros(s, n, m);
for t = 1 : s
    H  =  H0.* gamma(t);
    for l = 1 : n
        for k = 1 : m
            functionValue = functionQuadratic(X_SGD(:, 1), H);
            i = 1;
            while functionValue > functionLevel
                i = i + 1;
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * randn(d, 1); % Gaussian noise
%               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * 2 * sqrt(3) * (rand(d, 1)-0.5); % uniform noise
               X_SGD(:, i) = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sign(randn(d, 1)); % Bernoulli noise
%               noiseSph = randn(d, 1);
%               X_SGD(:, i)  = X_SGD(:, i - 1) - eta(l) * H * X_SGD(:, i - 1) + eta(l) * sigma * sqrt(d) * noiseSph/sqrt(noiseSph' * noiseSph); % uniform-spherical noise
               functionValue = functionQuadratic(X_SGD(:, i), H);
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
% save quadraticGauDat;
% save quadraticUniDat;
% save quadraticBerDat;
 save quadraticSphDat;

%% bugs and tips
% 1. any bug is possible - be brave and patient to debug

