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
%load data/ExampleData3
load

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


%% Maximum Likelihood (pooled)

MLEout = MLestimation(data,theta0,opts);

save MLEout

%% Plot Densities from any hierarchical parameters
if any(opts.Hier)
   MLEout.parh(opts.Hier) 
   
   d=gampdf(0:.01:1, MLEout.parh(opts.Hier==1),MLEout.parh(end-sum(opts.Hier)+1:end));
   %d=logncdf(0:.01:4, MLEout.parh(opts.Hier==1),MLEout.parh(end-sum(opts.Hier)+1:end));
   figure(1)
   plot(0:.01:1,d)
end
