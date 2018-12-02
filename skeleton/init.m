% function [mu,sigma,R,Q,Lambda_M] = init()
% This function initializes the parameters of the filter.
% Outputs:
%			mu(0):			3X1
%			sigma(0):		3X3
%			R:				3X3
%			Q:				2X2
function [mu,sigma,R,Q,Lambda_M] = init()
mu = [0;0;0]; % initial estimate of state
sigma = 1e-10*diag([1 1 1]); % initial covariance matrix
delta_m = 0.999;
Lambda_M = chi2inv(delta_m,2);

R = diag([0.01 0.01 0.0175]).^2;
Q = diag([0.01 0.0175]).^2;

%%% mean error(x, y, theta)=(0.000237, 0.000209, 0.000498)
%%% mean absolute error=(0.003173, 0.003925, 0.002562)
%%% total_time =18.852478

% R = diag([0.01 0.01 0.0175]).^2;
% Q = diag([0.2 0.2]).^2;

% mean error(x, y, theta)=(0.007140, 0.004798, -0.012870)
% mean absolute error=(0.041982, 0.042101, 0.044488)
% total_time =47.038265

% R = diag([1 1 1]).^2;
% Q = diag([0.1 0.1]).^2;
% mean error(x, y, theta)=(5.835831, 14.319903, -0.876833)
% mean absolute error=(16.127470, 17.422363, 1.417805)
% total_time =3.798329

% Using batch update
% mean error(x, y, theta)=(-0.026019, -0.022660, 0.018390)
% mean absolute error=(0.080351, 0.088318, 0.047883)
% total_time =6.649791


end