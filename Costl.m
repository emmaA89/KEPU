% -------------------------------------------------------------------------------- %
%      This software is licensed by Creative Commons BY-NC-SA:    %
%  http://creativecommons.org/licenses/by-nc-sa/3.0/it/legalcode %
% -------------------------------------------------------------------------------- %
%
% File: Costl(l,rbf,DM_data,Ytrain,noise,snr_data)
%
% Goal: script that defines the objective function for finding the 
%           optimal length scale paramter
%
% Inputs:   l:               variable to minimize
%              rbf:            radial basis function
%              DM_data:  training distance matrix
%              Ytrain:       the function values
%              noise:        noise in the measurments (either noise_F, i.e
%                                no noise, or noise_T)
%              snr_data:   signal to noise ratio
%
% Outputs:  cl:           the optimal value of the lenght scale
%
% -------------------------------------------------------------------------------- %
function cl = Costl(l,rbf,DM_data,Ytrain,noise,snr_data)
% Initialize and compute the kernel matrix
N = size(DM_data,1); K = rbf(l,DM_data); 
if noise == "noise_T"
    K = K + (1/snr_data)*eye(N);
end
% Compute the objective function
logdetK = (log(abs(det(K))));
c = K\Ytrain;
hsnorm = (Ytrain'*c);
cl = N*log(hsnorm) + logdetK;
end