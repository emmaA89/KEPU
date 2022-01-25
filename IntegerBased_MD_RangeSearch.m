%-------------------------------------------------------------------------%
%         This software is licensed by Creative Commons BY-NC-SA:         %
%      http://creativecommons.org/licenses/by-nc-sa/3.0/it/legalcode      %
%-------------------------------------------------------------------------%
%
% File: IntegerBased_MD_RangeSearch(puctr,puradius,dsites,index)
%
% Goal: find the data sites located in a given subdomain and the distances
%       between the subdomain centre and data sites
%
% Inputs:  puctr:       subdomain centre
%             puradius:  radius of PU subdomains
%             dsites:      NXM matrix representing a set of N data sites
%             index:       vector containing the indices of the data points 
%                             located in the k-th block (the block containing the 
%                             subdomain centre) and in the neighbouring blocks
%  
% Outputs: idx:  vector containing the indices of the data points located
%                       in a given PU subdomain
%               dist: vector containing the distances between the data sites 
%                       and the subdomain centre
%
%-------------------------------------------------------------------------%
function [idx, dist] = IntegerBased_MD_RangeSearch(puctr,puradius,...
    dsites,index)
N = size(dsites,1); dist = []; idx = []; % Initialize
% Compute distances between the data sites and the centre
for i = 1:N
    dist1(i) = norm(puctr - dsites(i,:));
end
% Use a sort procedure to order distances
if N > 0
    [sort_dist,IX] = sort(dist1);
    N1 = size(sort_dist,2); j1 = 1; j2 = 1; %Initialize
    % Find the data sites located in the given subdomain
    if nargin == 3
        idx = IX; dist = dist1;
    else
        while (j2 <= N1) && (sort_dist(j2) <= puradius)
            idx(j1) = index(IX(j2)); dist(j1) = dist1(IX(j2));
            j1 = j1 + 1; j2 = j2 + 1;
        end
    end
end