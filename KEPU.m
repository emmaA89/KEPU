% -------------------------------------------------------------------------------- %
%      This software is licensed by Creative Commons BY-NC-SA:    %
%  http://creativecommons.org/licenses/by-nc-sa/3.0/it/legalcode %
% -------------------------------------------------------------------------------- %
%
% File: KEPU(M,Xtrain,neval,npu,...
%            rbf,wf,Ytrain,method,kernel,noise,optim,l,sigma)
%
% Goal: script that performs partition of unity Kriging estimation
%
% Inputs:   M:             space dimension
%              Xtrain:        NXM matrix representing a set of N training data
%              neval:        number of test data in one direction
%              npu:           number of PU subdomains in one direction
%              rbf:             radial basis function
%              wf:              weight function
%              Ytrain:        the function values
%              method:     if set to 'canonical' it performs a canonical
%                                 implementation, else it uses the fitrgp.m routine
%              kernel:        keyword for the fitrgp.m function
%              noise:         noise in the measurments (either noise_F, i.e
%                                 no noise, or noise_T)
%              optim:        if set to 'opt_T', the KEPU optimizes the kernel
%                                parameters else it uses default parameters
%                                (opt_F) or user parameters (opt_U)
%              l:                length-scale parameter (optional)
%              sigma:       process variance (optional)
%
% Outputs:  KEPU_fit:           the KEPU fit
%                 KEPU_var:         the PU Kriging vairance
%                 Xtest:                test data
%
% Calls on: IntegerBased_MD_Structure,
%                IntegerBased_MD_Neighbourhood,
%                IntegerBased_MD_RangeSearch,
%                IntegerBased_MD_ContainingQuery,
%                DistanceMatrix, MakeSDGrid
%
% Remark:   DistanceMatrix and MakeSDGrid come from the book:
%           [G. E. Fasshauer, Meshfree approximation methods with
%           Matlab, World Scientific, Singapore, 2007].
%                 IntegerBased_MD_Structure,
%                 IntegerBased_MD_Neighbourhood,
%                 IntegerBased_MD_RangeSearch,
%                 IntegerBased_MD_ContainingQuery, are  by
%           R. Cavoretto, A. De Rossi, E. Perracchione, are avilable at
%           available at: http://hdl.handle.net/2318/1559094
%
% Required toolbox: Signal processing, optimization and
%                              Statistics and Machine Learning Toolbox
%
% -------------------------------------------------------------------------------- %
function [KEPU_fit, KEPU_var, Xtest] = KEPU(M,Xtrain,neval,npu,...
    rbf,wf,Ytrain,method,kernel,noise,optim,l,sigma)
% Check the input parameters for the kernel.
if optim == "opt_U" && (isempty(l) || isempty(sigma))
    disp('Error: cannot specify optim_U if sigma and l are empty')
    return
end
if  optim == "opt_T" && noise == "noise_T"
    % Compute the signal to noise ratio (used
    % to smooth out the noise)
    snr_data = abs(snr(Ytrain));
end
% Create npu^M equally spaced PU centres
puctrs = MakeSDGrid(M,npu);
% Create neval^M equally spaced test data
Xtest = MakeSDGrid(M,neval);
puradius = sqrt(2)/npu;  % Define the initial PU radius
wep = 1/puradius;  % Parameter for weight function
npu_M = size(puctrs,1); neval_M = size(Xtest,1); % Initialize
KEPU_fit = zeros(neval_M,1);  % Initialize
KEPU_var = zeros(neval_M,1);  % Initialize
% Compute the Shepard evaluation matrix and its square
DM_eval = DistanceMatrix(Xtest,puctrs);
SEM = wf(wep,DM_eval);
SEM = spdiags(1./(SEM*ones(npu_M,1)),0,neval_M,neval_M)*SEM;
SEM_sigma = spdiags(1./(SEM.^2*ones(npu_M,1)),0,neval_M,...
    neval_M)*SEM.^2;
% Parameter for the integer-based partitioning structure
q = ceil(1./puradius);
% Build the partitioning structure for training and test data
idx_train = IntegerBased_MD_Structure(Xtrain,q,puradius,M);
idx_test = IntegerBased_MD_Structure(Xtest,q,puradius,M);
if optim == "opt_F" % Define default parameters
    k = 1; % Initialize
    for j = 1:npu_M % Loop over subdomains
        % Find the block containing the j-th subdomain centre
        index1{j} = IntegerBased_MD_ContainingQuery(puctrs(j,:),q,...
            puradius,M);
        % Find the training data located in the j-th subdomain
        [dxx, dx] = IntegerBased_MD_Neighbourhood(Xtrain,idx_train,...
            index1{j},q,M,1);
        idx1{j} = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,...
            dxx,dx);
        if ~isempty(idx1{j})  % Compute the local default parameters
            lvec(k) = mean(std(Xtrain(idx1{j},:)));
            sigmavec(k) = std(Ytrain(idx1{j}))/sqrt(2);
            k = k+1;
        end
    end
    % Define the default parameters
    l = mean(lvec);
    sigma = mean(sigmavec);
end
% Main loop for KEPU
for j = 1:npu_M % Loop over subdomains
    if optim == "opt_F"
        % Define the training data located in the j-th subdomain
        idx = idx1{j};
        % Find the test data located in the j-th subdomain
        [edxx, edx] = IntegerBased_MD_Neighbourhood(Xtest,idx_test,...
            index1{j},q,M,1);
        eidx = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,...
            edxx,edx);
    else
        % Find the block containing the j-th subdomain centre
        index = IntegerBased_MD_ContainingQuery(puctrs(j,:),q,...
            puradius,M);
        % Find the training data located in the j-th subdomain
        [dxx, dx] = IntegerBased_MD_Neighbourhood(Xtrain,idx_train,...
            index,q,M,1);
        idx = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,...
            dxx,dx);
        % Find the test data located in the j-th subdomain
        [edxx, edx] = IntegerBased_MD_Neighbourhood(Xtest,idx_test,...
            index,q,...
            M,1);
        eidx = IntegerBased_MD_RangeSearch(puctrs(j,:),puradius,...
            edxx,edx);
    end
    if  (~isempty(eidx)) &&  (~isempty(idx))
        if  optim == "opt_T" || method == "canonical"
            % Compute the local distance matrix
            DM_data = DistanceMatrix(Xtrain(idx,:),Xtrain(idx,:));
            if  optim == "opt_T"
                % Define a range for the length scale parameter
                min_l = mean(std(Xtrain(idx)))/10;
                max_l = mean(std(Xtrain(idx)))*10;
                if noise == "noise_T"
                    % Find the optimal  length scale parameter
                    l =  fminbnd(@(l)Costl(l,rbf,DM_data,Ytrain(idx),...
                        noise,snr_data),min_l,max_l);
                else
                    l =  fminbnd(@(l)Costl(l,rbf,DM_data,Ytrain(idx),...
                        noise,0),min_l,max_l);
                end
            end
            K = rbf(l,DM_data); % Interpolation matrix
            if noise == "noise_F"
                % Compute the interpolation coefficients
                coef = K\Ytrain(idx);
            else
                % Compute the interpolation coefficients via
                % regularization
                coef =  (K + 1/snr_data*eye(size(K)))\Ytrain(idx);
            end
        end
        if optim == "opt_T" % Define the optimal process variance
            sigma = sqrt(1/length(idx)*Ytrain(idx)'*coef);
        end
        if method == "fitrgp" % Call fitrgp.m and fit the local model
            Mdl = fitrgp(Xtrain(idx,:),Ytrain(idx),'KernelParameters',...
                [l,sigma],'KernelFunction',kernel);
            % Evaluate the local model and define the
            % local variance
            [localfit, std_fit] = predict(Mdl,Xtest(eidx,:));
            loc_var =  std_fit.^2/4;
            if noise == "noise_T"
                % Compute the fit on the training data (used to
                % define the confidence intervals)
                localfit_data = predict(Mdl,Xtrain(idx,:));
            end
        else
            % Compute the evaluation matrix and the local fit
            DM_eval = DistanceMatrix(Xtest(eidx,:),Xtrain(idx,:));
            K_eval = rbf (l, DM_eval); localfit = K_eval*coef;
            if noise == "noise_T"
                % Compute the fit on the training data (used to
                % define the confidence intervals)
                localfit_data = K*coef;
            end
            % Compute the local variance
            loc_var =  sigma^2*(rbf(l,0) - sum((K_eval/K).*K_eval,2));
        end
        % Accumulate the global fit and variance
        KEPU_fit(eidx) = KEPU_fit(eidx) + localfit.*SEM(eidx,j);
        if noise == "noise_F"
            KEPU_var(eidx) = KEPU_var(eidx) +  loc_var.*...
                SEM_sigma(eidx,j);
        else
            % Add the estimated noise to the Kriging variance
            noise_est = max(abs(localfit_data-Ytrain(idx)));
            KEPU_var(eidx) = KEPU_var(eidx) +  (loc_var+noise_est^2/4)....
                .*SEM_sigma(eidx,j);
        end
    end
end