% parameters
%% R0 um pouco menor do que 1 %%
% parameters
N=100;
R0 = 0.95;
mu = 1/(75*365); %mortality
gamma = 1/10; %recovery rate
beta = R0*(mu+gamma)/N; %infection rate
phi = 0.8; %vaccination rate
sigma = 0.25; % vaccination 'failure'
theta = 1/(365*3); % wanning vaccine immunity
psi = 1/180; % wanning natural imunity
parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi);
%%
% another set of parameters
N = 100;
mu = 0.00004; %mortality
gamma = 1/7; %recovery rate
beta = R0*(mu+gamma)/N; %infection rate
%beta = 1/6;
phi = 0.9; %vaccination rate
sigma = 0.05; % vaccination 'failure'
theta = 0.9; % wanning vaccine immunity
psi = 1/180; % wanning natural imunity
parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi);
%% Abdu Parameters
%mu = 0.0000365, gamma = 1/7, psi = 1/90, p = 0.9, theta = 1/1095, varepsilon = 0.8, sigma = 1 - 0.8, R0 = 3.5, N = 100
N = 100;
mu = 0.0000365; %mortality
gamma = 1/7; %recovery rate
R0 = 3.5;
beta = R0*(mu+gamma)/N; %infection rate
theta = 1/1095; % wanning vaccine immunity
p = 0.90;
phi = p*(mu+theta)/(1-p); %vaccination rate
%phi = 0;
sigma = 1 - 0.8; % vaccination 'failure'
psi = 1/90; % wanning natural imunity
parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi);
%% Modified Abdu Parameters
%mu = 0.0000365, gamma = 1/7, psi = 1/90, p = 0.9, theta = 1/1095, varepsilon = 0.8, sigma = 1 - 0.8, R0 = 3.5, N = 100
N = 100;
mu = 0.0000365; %mortality
gamma = 1/7; %recovery rate
R0 = 3.5;
beta = R0*(mu+gamma)/N; %infection rate
theta = 1/1095; % wanning vaccine immunity
p = 0.913; %p = 0.893;
phi = p*(mu+theta)/(1-p); %vaccination rate
%phi = 0;
sigma = 1 - 0.8; % vaccination 'failure'
psi = 1/90; % wanning natural imunity
parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi);
%% Another set of Abdul Parameters
%mu = 0.0000365, gamma = 1/7, psi = 1/90, p = 0.9, theta = 1/1095, varepsilon = 0.8, sigma = 1 - 0.8, R0 = 3.5, N = 100
N = 100;
mu = 0.0000365; %mortality
gamma = 1/7; %recovery rate
R0 = 3.5;
beta = R0*(mu+gamma)/N; %infection rate
theta = 1/1095; % wanning vaccine immunity
psi = 1/90; % wanning natural imunity
p = (mu+psi)*(R0-1)/((mu+psi)*R0+theta-psi)+0.001;
p = 25/28;
p = 0.91269625;
phi = p*(mu+theta)/(1-p); %vaccination rate
%phi = 0;
sigma = 1 - 0.8; % vaccination 'failure'
parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi);
%%
function parameters = SetParameters(N, mu, gamma, beta, theta, phi, sigma, psi)

    parameters = [];
    parameters.N = N;
    parameters.mu = mu; %mortality
    parameters.gamma = gamma; %recovery rate
    parameters.beta = beta; %infection rate
    parameters.theta = theta; % wanning vaccine immunity
    parameters.phi = phi; %vaccination rate
    parameters.sigma = sigma; % vaccination 'failure'
    parameters.psi = psi; % wanning natural imunity
end