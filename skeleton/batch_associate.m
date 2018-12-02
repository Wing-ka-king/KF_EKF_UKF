% function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           Q                   2X2
%           z(t)                2Xn
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1Xn
%           outlier             1Xn
%           nu_bar(t)           2nX1
%           H_bar(t)            2nX3
function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,M,Lambda_m,Q)
    N1 = size(z);
    N2 = size(M);
    nu_bar = zeros(2,1);
    H_bar = zeros(2,3);
    for i = 1:N1(2)
        for j = 1:N2(2)
                z_hat = observation_model(mu_bar,M,j);
                H(:,:,j) = jacobian_observation_model(mu_bar, M, j, z_hat, 1);
                S(:,:,j) = H(:,:,j)*sigma_bar*H(:,:,j).' + Q;
                nu(:,j) = z(:,i) - z_hat;
                nu(2,j) = mod(nu(2,j)+pi,2*pi)-pi;
                D(:,j) = nu(:,j).'*inv(S(:,:,j))*nu(:,j);
                psi(j) = (det(2*pi*S(:,:,j)))^(-1/2)*exp((-1/2)*D(:,j));
        end
        
        c(i) = find(psi==max(psi));
        outlier(i) = (D(:,c(i)) >= Lambda_m);
        nu_bar(:,i) = nu(:,c(i));    
        H_bar_temp(:,:,i) = H(:,:,c(i));
    
    end
    nu_bar = reshape(nu_bar, 2*N1(2),1);
    H_bar  = H_bar_temp(:,:,1);
    for i = 2:N1(2)
        H_bar = [H_bar;H_bar_temp(:,:,i)];
    end

end

