%% Network regularised Cox Regession
clc
clear all

load Ovarian

Data = Normal(Data);

S = getNormalizedNetworkInformation(Data);

lambda = [1e-5,1e-4];

alpha = [0.01,0.5,0.95];
beta_coexpression = NetworkRegularisedCox(Data,lambda,alpha,d, S);

