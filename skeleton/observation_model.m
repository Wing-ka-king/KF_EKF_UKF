% function h = observation_model(x,M,j)
% This function is the implementation of the h function.
% The bearing should lie in the interval [-pi,pi)
% Inputs:
%           x(t)        3X1
%           M           2XN
%           j           1X1
% Outputs:  
%           h           2X1
function h = observation_model(x,M,j)
h(1,1) = sqrt((M(1,j)-x(1))^2+(M(2,j)-x(2))^2);
h(2,1) = atan2(M(2,j)-x(2), M(1,j)-x(1))-x(3);
h(2,1)= mod(h(2,1)+pi,2*pi)-pi;
end