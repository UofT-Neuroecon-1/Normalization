%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation / Model-comparison Example File %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

addpath(['..' filesep 'EstimationAdaptive' filesep])
addpath(['DERIVESTsuite' filesep])

%% Load Backup
% if you want to load a backup from a previous estimation (set to [] if
% not)
backup_file = ''; 'backup1-StndPDN.mat';

%% Data format
% Requires a struct array (e.g. data) of size N x 1
% Each entry must contain a structure with at least 2 fields:
% data(n).X : cell array of T observations.
% data(n).X{t}: Matrix of J(t) options x K(t) attributes
%    -assumes that max(J(t)) is the same across all subjects (i.e. each
%    subject sees the biggest choice set)
% data(n).y : Vector of T x 1 Choices

%load data/ExampleDataSS
load data/ExampleData3

%% Estimation Parameters
% You can create your own parameters that will be passed to the likelihood
% function and the particles initialization functions

% Models to use for estimations
% Each model should have a corresponding entry in the following files:
% InitParticle : Returns a draw from prior for one particle
% Mutate : The Metropolis-Hastings Mutation step
% ProbaChoice : The likelihood of one observation given a model and particle

%{'MNP'}: MNP
%{'DNw'}: allow asymmetric weights
%{'DNwb'}: allow asymmetric weights and beta
%{'DN'}: basic DN
%{'DNb'}: basic DN with beta
%{'PDNNew'}
%{'RemiStand';'HierarchicalProbit'}
%{'Range'}
opts.Models = {'DNb'}

opts.Prob='Probit'; %'Probit','Logit','GHK', 'HP' Is the covariance matrix restricted to be independent?
opts.names={'kappa','sigma','omega','a','b','w2'};
opts.i0=3; %normalize w.r.t. altnerative...
opts.Hier=logical([0 0 0 0 0 0]); %which parameters to make hierarchical: kappa, sigma, omega, alpha, beta;

%opts.cluster=ones(length(data),1); %all in same cluster (i.e. pooled)
opts.cluster=length(data); % each in own cluster

% opts.Sinz=0;
% opts.weights=0;

%%%%Initial Parameter Values (the number of elements must match the number
%%%%of parameters needed in the model above, +1 for every hierarchical
%%%%parameter)
%theta0=[0.027 0.214 3.458 0.079];
%theta0=[0.012 0.412 513.7 0];
%theta0=[0.012 0.412 25.74];
%theta0=[0 0.5 51 0];
%theta0=[0.0265 0.214 3.4582 0.0792];
%theta0=[];
%theta0=[0 0.4923 19 -.032];
%theta0=[0.012 0.412 25.74];
%theta0=[0.012 0.412];
theta0=[];
%theta0=[0.069281675331641 1.823034213307364 0.108391494341653]; %random parameter on omega
%theta0=[0.069281675331641 1.823034213307364 1 0.108391494341653]; %random parameter on omega w beta
%theta0=[0.441 0.357];
%theta0=[0.28346 0.61467 1.1277]; %set size: Range. random parameter on omega


%%%%Estimation Specific Parameters
%for particle filter
opts.G = 3; % Number of particles group
opts.P = 128; % Number of particles per group
opts.Adaptive = true; % Use the adaptive SMC (see Durham, Geweke 2014).
opts.Msteps = 10; % Number of mutate steps
opts.ress_threshold = 0.8;

%for ML
opts.getP=0; %Set this to 1 if you just want to evaluate at initial parameter vector
opts.numInit=1; %1: run using gradient- method first at theta0. If >1, use random starting points with gradient free method.
opts.R=5000;
opts.ses=1; %calculate ses? Set to zero if you just want estimates (much faster)
opts.Display='iter-detailed';

%%%% Model specific parameters
opts.NormDraw = mvnrnd(zeros(4,1),eye(4),1000); % Pre-draw errors for GHK (if needed)

opts.attrVals{1}= (0:4);          opts.attrNames{1}="Bid";
opts.attrSign = [1];

opts.K = numel(opts.attrVals); 
opts.attrMax = zeros(1,opts.K);
for k=1:opts.K
    opts.attrMax(k) = max(opts.attrVals{k});
end

clear k

%Save format
opts.Tag = 'StndPDN'; % This tag will be added to the output file
opts.savefile = ['Analysis' filesep opts.Tag sprintf('-%.0fx%.0f-M%.0f-',opts.G,opts.P,opts.Msteps) datestr(datetime('now'),'yyyy-mm-dd-HH.MM') '.mat'];
if ~exist('Analysis','dir')
    mkdir('Analysis')
end


 %% Estimation: use adaptive algorthm

% EstimationOutput = EstimationAdaptiveSMC( datapooled, opts, backup_file);
% Particles = EstimationOutput.Particles;
% 
% if strcmp(opts.Models{1},'PDNNew')
%     EstimationOutput.Particles{1}.postmeans
%     fprintf('Posterior Mean (Across Subjects) \n')
%     fprintf('alpha: %f \n sigma: %f \n omega: %f \n',mean(EstimationOutput.Particles{1}.postmeans))
% end

% sample_part = EstimationOutput.Particles(1).particle{1,1};
% vect_theta = zeros(opts.G,opts.P,size(sample_part.theta,1),size(sample_part.theta,2));
% for g=1:opts.G
%     for p=1:opts.P
%         vect_theta(g,p,:,:) = EstimationOutput.Particles(1).particle{g,p}.theta;
%     end
% end
% subplot(2,1,1);
% histogram(vect_theta(:,:,1,1),0:0.05:2)
% subplot(2,1,2);
% histogram(vect_theta(:,:,1,2),0:0.05:2)
 %% Full posterior (all the subjects superposed)
% prior = [betarnd(3,1,10000,1) gamrnd(1,0.5,10000,1) gamrnd(1,1,10000,1)];
% prior_mean = mean(prior,1);
% % Compute posterior means
% size_NK = size(Particles{1}.particle{1}.theta);
% VectorizedTheta = nan(opts.P*opts.G,size_NK(1),size_NK(2));
% for p=1:opts.G*opts.P
%     VectorizedTheta(p,:,:) = Particles{1}.particle{p}.theta;
% end
% post_mean = squeeze(mean(mean(VectorizedTheta,1),2,'omitnan'))
% % Plot posteriors
% for theta=1:3
%     subplot(2,3,theta)
%     histogram(prior(:,theta),20,'Normalization','pdf')
%     title(sprintf('prior mean: %.2f',prior_mean(theta)));
%     subplot(2,3,3+theta)
%     hist_data = VectorizedTheta(:,:,theta);
%     hist_data = hist_data(~isnan(hist_data(:)))
%     histogram(hist_data(:),20,'Normalization','pdf')
%     title(sprintf('post mean: %.2f',post_mean(theta)));
% end

  %% Maximum Likelihood (within)
% theta0 = [1,0,.2];
% options = optimoptions('fminunc','OptimalityTolerance',1.0e-9,'Display','iter-detailed');
% theta_indiv_ml = nan(numel(SubjData),3);
% for subj = 1:numel(SubjData)
%     Target = @(theta) -LogLikelihood( SubjData{subj}.Xs, SubjData{subj}.ChoiceList, 1 , 'DN' , theta, opts );
%     theta_indiv_ml(subj,:) = fminunc(Target,theta0,options);
% end
%% Maximum Likelihood (pooled)
% pool data


% theta0 = [1,0,.2];
% options = optimoptions('fminunc','OptimalityTolerance',1.0e-9,'Display','iter-detailed');
% Target = @(theta) -LogLikelihood( PooledXs, PooledChoiceList, 1 , 'DN' , theta, opts );
% 
% theta_pooled = fminunc(Target,theta0,options)

% data_s1 = struct;
% data_s1.X = datapooled.X(1:270);
% data_s1.y = datapooled.y(1:270);
% data_s1.J = datapooled.J(1:270);
% data_s1.K = datapooled.K;
% MLEout = MLestimation(data_s1,[],opts);

MLEout = MLestimation(data,theta0,opts);

save MLEout
%% Plot Likelihoods conditional on true values
% true_theta = [par.a,par.sigma,par.omega]
% gridpoints = linspace(0,3);
% yy = nan(numel(gridpoints),3);
% for th = 1:3
%     theta = true_theta;
%     for i = 1:numel(gridpoints)
%         theta(th) = gridpoints(i);
% 
%         yy(i,th) = LogLikelihood( PooledXs, PooledChoiceList, 1 , 'DN' , theta, opts );
% 
%     end
%     subplot(3,1,th);
%     plot(gridpoints,yy(:,th));
%     title(sprintf('LogLik param %d',th));
% end

%%
if any(opts.Hier)
   MLEout.parh(opts.Hier) 
   
   d=gampdf(0:.01:1, MLEout.parh(opts.Hier==1),MLEout.parh(end-sum(opts.Hier)+1:end));
   %d=logncdf(0:.01:4, MLEout.parh(opts.Hier==1),MLEout.parh(end-sum(opts.Hier)+1:end));
   figure(1)
   plot(0:.01:1,d)
end
