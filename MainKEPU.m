% -------------------------------------------------------------------------------- %
%      This software is licensed by Creative Commons BY-NC-SA:    %
%  http://creativecommons.org/licenses/by-nc-sa/3.0/it/legalcode %
% -------------------------------------------------------------------------------- %
%
% The following script provides the 1d KEPU examples included in
% R. Cavoretto, A. De Rossi, E. Perracchione, 
% Learning with Partition of Unity-based Kriging Estimators: 
% Efficient Uncertainty Computation. Available on Researchgate at
% https://www.researchgate.net/profile/Emma-Perracchione
%
% -------------------------------------------------------------------------------- %

close all
clear all
warning off

% NOTE: to save results, create in the current Matlab path a
% folder named Results

% ----------------- The following functions are needed: ----------------- %
% DistanceMatrix and MakeSDGrid, from the book
%        [G. E. Fasshauer, Meshfree approximation methods with
%        Matlab, World Scientific, Singapore, 2007].
%  IntegerBased_MD_Structure, IntegerBased_MD_Neighbourhood,
%        IntegerBased_MD_RangeSearch, IntegerBased_MD_ContainingQuery,
%        by R. Cavoretto, A. De Rossi, E. Perracchione, avilable at
%        available at: http://hdl.handle.net/2318/1559094
% shade, by Javier Montalt Tordera, available at
%       https://it.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot
% haltonseq, by D. Dougherty, available at
%       http://www.math.usm.edu/cschen/COS702/HW1/Reese/haltonseq.m
%
% Required toolbox: Signal processing, optimization and
%                              Statistics and Machine Learning Toolbox

% ----------------------------------------------------------------------- %
Krig = "KEPU"; % Keyword for saving results
% Define the number of training data
Nvec = [9,9,17,33,65,129].^2;
% IMPORTANT NOTE: the first test is repeated twice to avoid
% misleading CPU times due to the opening of the
% graphics tools
% Define the number of patches in one direction
npuvec = floor(sqrt(Nvec));
% Initialize
CPUTime = zeros(1,length(Nvec)-1); % KEPU CPU times
RMSE = zeros(1,length(Nvec)-1); % KEPU RMSE
AKV = zeros(1,length(Nvec)-1); % KEPU Average Kriging
% Variance
MAE = zeros(1,length(Nvec)-1); % KEPU Maximum Absolute
% Error
MKV = zeros(1,length(Nvec)-1); % KEPU Maximum Kriging
% Variance
M = 1; % Problem dimension
Ms = num2str(M); % Keyword for saving results
neval = 100; % Number of test data in one direction

% -------    Uncomment one of the following 4 tests  -------- %

% -----------------------------  Example 1 -------------------------- %
%
% type = "Halton"; % Data type. In this script it can be
% % either Halton, Equispaced (that leads to grids in the
% % multivariate case) or Random.
% kernel = "matern32";  % Kernel type. In this script it
% % can be either either matern32
% % (Matern C^2),  matern52 (Matern C^4) or
% % squaredexponential (the Gaussian)
% method = "canonical";  % Implementation type (either,
% % canonical or fitrgp)
% optim ="opt_F"; % Parameter for setting the kernels length
% % scale (l) and process variance (sigma). It can be either
% % opt_F (default parameters), opt_T (optimized parameters)
% % or opt_U (user parameters).
% % if optim  is set to opt_U then specify l and sigma
% % else do as follows
% if optim == "opt_F" || optim == "opt_T"
%     l =[]; sigma = [];
% end
% noise ="noise_F"; % Noise (either noise_F (no noise) or
% % noise_T (add noise to measurments))
% f_t = @(x) sin(10*pi.*(x+0.1)); % Define the test function
% a = 0; b = 1; % Omega
% % Note: Data should be rescaled in [0,1]^d

% % -----------------------------  Example 2 -------------------------- %
%
% type = "Halton"; % Data type. In this script it can be
% % either Halton, Equispaced (that leads to grids in the
% % multivariate case) or Random.
% kernel = "matern32";  % Kernel type. In this script it
% % can be either either matern32
% % (Matern C^2),  matern52 (Matern C^4) or
% % squaredexponential (the Gaussian)
% method = "canonical";  % Implementation type (either,
% % canonical or fitrgp)
% optim ="opt_T"; % Parameter for setting the kernels length
% % scale (l) and process variance (sigma). It can be either
% % opt_F (default parameters), opt_T (optimized parameters)
% % or opt_U (user parameters).
% % if optim  is set to opt_U then specify l and sigma
% % else do as follows
% if optim == "opt_F" || optim == "opt_T"
%     l =[]; sigma = [];
% end
% noise ="noise_F"; % Noise (either noise_F (no noise) or
% % noise_T (add noise to measurments))
% f_t = @(x) sin(10*pi.*(x+0.1)); % Define the test function
% a = 0; b = 1; % Omega
% % Note: Data should be rescaled in [0,1]^d

%
% % -----------------------------  Example 3 -------------------------- %
%
% type = "Equispaced"; % Data type. In this script it can be
% % either Halton, Equispaced (that leads to grids in the
% % multivariate case) or Random.
% kernel = "squaredexponential";  % Kernel type. In this script it
% % can be either either matern32
% % (Matern C^2),  matern52 (Matern C^4) or
% % squaredexponential (the Gaussian)
% method = "fitrgp";  % Implementation type (either,
% % canonical or fitrgp)
% optim ="opt_F"; % Parameter for setting the kernels length
% % scale (l) and process variance (sigma). It can be either
% % opt_F (default parameters), opt_T (optimized parameters)
% % or opt_U (user parameters).
% % if optim  is set to opt_U then specify l and sigma
% % else do as follows
% if optim == "opt_F" || optim == "opt_T"
%     l =[]; sigma = [];
% end
% noise ="noise_T"; % Noise (either noise_F (no noise) or
% % noise_T (add noise to measurments))
% noise_level = 0.1; % Parameter for Gaussian white noise
% % Define the test function
% f_t = @(x) cos(14*pi.*(x+0.5)) ./ (2.*x+0.5)+(x+0.5-1).^4;
% a = 0; b = 1; % Omega
% % Note: Data should be rescaled in [0,1]^d

%
% % -----------------------------  Example 4 -------------------------- %
type = "Random"; % Data type. In this script it can be
% either Halton, Equispaced (that leads to grids in the
% multivariate case) or Random.
kernel = "matern52";  % Kernel type. In this script it
% can be either either matern32
% (Matern C^2),  matern52 (Matern C^4) or
% squaredexponential (the Gaussian)
method = "fitrgp";  % Implementation type (either,
% canonical or fitrgp)
optim ="opt_T"; % Parameter for setting the kernels length
% scale (l) and process variance (sigma). It can be either
% opt_F (default parameters), opt_T (optimized parameters)
% or opt_U (user parameters).
% if optim  is set to opt_U then specify l and sigma
% else do as follows
if optim == "opt_F" || optim == "opt_T"
    l =[]; sigma = [];
end
noise ="noise_T"; % Noise (either noise_F (no noise) or
% noise_T (add noise to measurments))
noise_level = 0.1; % Parameter for Gaussian white noise
% Define the test function
f_t = @(x) cos(14*pi.*(x+0.5)) ./ (2.*x+0.5)+(x+0.5-1).^4;
a = 0; b = 1; % Omega
% Note: Data should be rescaled in [0,1]^d

%
% ------------- Define the KEPU setting ------------------ %

% Define the PU weights (Wendland C^2)
wf = @(e,r) (max(1-(e*r),0).^4).* ...
    (4*(e*r)+1);

% Define the kernel
if kernel == "matern32"
    rbf = @(l,r) (1+sqrt(3)/l.*r).*exp(-sqrt(3)/l.*r);
else
    if kernel == "squaredexponential"
        rbf = @(l,r) exp(-r.^2/(2.*l.^2));
    else
        rbf = @(l,r) (1+sqrt(5)/l.*r+(5*r.^2)./(3*l))...
            .*exp(-sqrt(5)/l.*r);
    end
end

% ----------------------------  Loop over n  ------------------------------ %
for k = 1:length(Nvec)
    N = Nvec(k); % Number of data
    % Define the training data
    if type == "Halton"
        Xtrain = haltonseq(N,M);
    else
        if type == "Equispaced"
            N_1 = floor(N^(1/M));
            Xtrain = MakeSDGrid(M,N_1);
            N = size(Xtrain,1);
        else
            rng(10) % Fix the random seed
            Xtrain = rand(N,M);
        end
    end
    % Define the training samples
    if noise == "noise_F"
        Ytrain = f_t(Xtrain);
    else
        rng(10) % Fix the random seed and perturb data
        Ytrain = f_t(Xtrain)+noise_level*randn(N,1);
    end
    npu = npuvec(k); % Define the number of
    % PU patches
    
    tic % Call KEPU and compute CPU times
    [KEPU_fit, KEPU_var, Xtest] = KEPU(M,Xtrain,neval, ...
        npu,rbf,wf,Ytrain,method,kernel,noise,optim,l,sigma);
    t1 = toc;
    
    if k >1
        k_new = k-1;
        % Print details
        disp('---------------------------------------------------------------')
        formatSpec = 'KEPU with N:                            %d\n';
        fprintf(formatSpec,N) % Print
        Ns = num2str(N); % Keywords for saving results
        all_param = strcat(Krig,Ms,'_',type,'_',kernel,'_', ...
            method,'_',optim,'_',noise,'_',Ns);
        
        % Compute CPU times and accuracy indicators
        Ytest = f_t(Xtest);
        CPUTime(k_new) = t1;
        RMSE(k_new) = sqrt(norm((KEPU_fit-Ytest),2)^2/N);
        AKV(k_new) = mean(KEPU_var);
        MAE(k_new) = max(abs(KEPU_fit-Ytest));
        MKV(k_new) = max(KEPU_var);
        
        % Save and print reults
        Folder = "Results";
        UpFolder = pwd;
        currentFolderString = convertCharsToStrings(UpFolder);
        StringaOutput = fullfile(currentFolderString, Folder);
        PathOutput = convertStringsToChars(StringaOutput);
        
        if noise == "noise_F"
            formatSpec = 'CPUTime [s]:                               %f\n';
            fprintf(formatSpec,CPUTime(k_new))
            formatSpec = 'Root Mean Squared Error:           %e\n';
            fprintf(formatSpec,RMSE(k_new))
            formatSpec = 'Maximum Absolute Error:            %e\n';
            fprintf(formatSpec,MAE(k_new))
        else
            formatSpec = 'CPUTime [s]:                               %f\n';
            fprintf(formatSpec,CPUTime(k_new))
            formatSpec = 'Averaged of Kriging Variance :    %e\n';
            fprintf(formatSpec,AKV(k_new))
            formatSpec = 'Maximum Kriging Variance:         %e\n';
            fprintf(formatSpec,MKV(k_new))
        end
        
        % Print and save the figure for plotting
        h = figure(k_new);
        xlabel('x_1', 'FontSize',11)
        ylabel('f', 'FontSize',11)
        ax = gca;
        ax.FontSize = 11;
        axis square
        grid on
        box on
        hold on
        if noise == "noise_T"
            % Plot training data
            plot(Xtrain,Ytrain,'k.', 'MarkerSize',3)
        end
        % Plot the fit
        plot(Xtest,KEPU_fit,'-b', 'LineWidth',3)
        % Plot the confidence intervals
        shade(Xtest,KEPU_fit,Xtest,KEPU_fit+2*sqrt(KEPU_var),...
            'FillType',[1 2;2 1],'FillColor','m','FillAlpha',0.1);
        shade(Xtest,KEPU_fit,Xtest,KEPU_fit-2*sqrt(KEPU_var),...
            'FillType',[1 2;2 1],'FillColor','m','FillAlpha',0.1);
        all_param_old = all_param;
        all_param = strcat(all_param_old,'_','fit');
        tmp = fullfile(PathOutput,all_param);
        saveas(h,tmp,'png')
        saveas(h,tmp,'fig')
    end
end

% Save all variables
all_param_old = strcat(Krig,Ms,'_',type,'_',kernel,'_',...
    method,'_',optim,'_',noise);
temp=['CPUTime','.mat'];
all_param = strcat(all_param_old,'_',temp);
tmp = fullfile(PathOutput,all_param);
save(tmp,'CPUTime')
temp=['RMSE','.mat'];
all_param = strcat(all_param_old,'_',temp);
tmp = fullfile(PathOutput,all_param);
save(tmp,'RMSE')
temp=['AKV','.mat'];
all_param = strcat(all_param_old,'_',temp);
tmp = fullfile(PathOutput,all_param);
save(tmp,'AKV')
temp=['MAE','.mat'];
all_param = strcat(all_param_old,'_',temp);
tmp = fullfile(PathOutput,all_param);
save(tmp,'MAE')
temp=['MKV','.mat'];
all_param = strcat(all_param_old,'_',temp);
tmp = fullfile(PathOutput,all_param);
save(tmp,'MKV')

disp(' Summarizing results' )
if noise == "noise_F"
    disp('---------------------------------------------------------------')
    disp(' CPU times varying n ' )
    formatSpec = 'CPUTime [s]:                               %f\n';
    fprintf(formatSpec,CPUTime)
    disp('---------------------------------------------------------------')
    disp(' RMSE varying n ' )
    formatSpec = 'Root Mean Squared Error:           %e\n';
    fprintf(formatSpec,RMSE)
    disp('---------------------------------------------------------------')
    disp(' MAE varying n ' )
    formatSpec = 'Maximum Absolute Error:            %e\n';
    fprintf(formatSpec,MAE)
else
    disp('---------------------------------------------------------------')
    disp(' CPU times varying n ' )
    formatSpec = 'CPUTime [s]:                               %f\n';
    fprintf(formatSpec,CPUTime)
    disp('---------------------------------------------------------------')
    disp(' AKV varying n ' )
    formatSpec = 'Averaged of Kriging Variance :    %e\n';
    fprintf(formatSpec,AKV)
    disp('---------------------------------------------------------------')
    disp(' MKV varying n ' )
    formatSpec = 'Maximum Kriging Variance:         %e\n';
    fprintf(formatSpec,MKV)
end

% Save and print figures for sum up the results
h=figure(k_new+1);
all_param_figs = convertStringsToChars(...
    strcat(Krig,Ms,'_',type,'_',kernel,'_',...
    method,'_',optim,'_',noise,'_','CPUtimes'));
loglog(Nvec(2:end),CPUTime,'m.', 'MarkerSize',22)
xlabel('n', 'FontSize',11)
ylabel('CPU', 'FontSize',11)
ax = gca;
ax.FontSize = 11;
axis square
box on
grid on
hold on
loglog(Nvec(2:end),CPUTime,'m-', 'LineWidth',2)
temp=[all_param_figs,'.fig'];
saveas(h,fullfile(PathOutput, temp))
temp=[all_param_figs,'.epsc'];
saveas(h,fullfile(PathOutput, temp))
if noise == "noise_F"
    h=figure(k_new+2);
    all_param_figs = convertStringsToChars(...
        strcat(Krig,Ms,'_',type,'_',kernel,'_',...
        method,'_',optim,'_',noise,'_','RMSEs'));
    loglog(Nvec(2:end),RMSE,'m.', 'MarkerSize',22)
    xlabel('n', 'FontSize',11)
    ylabel('RMSE', 'FontSize',11)
    ax = gca;
    ax.FontSize = 11;
    axis square
    box on
    grid on
    hold on
    loglog(Nvec(2:end),RMSE,'m-', 'LineWidth',2)
    temp=[all_param_figs,'.fig'];
    saveas(h,fullfile(PathOutput, temp))
    temp=[all_param_figs,'.epsc'];
    saveas(h,fullfile(PathOutput, temp))
end

