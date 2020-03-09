%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation / Model-comparison Example File %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

addpath(['..' filesep 'EstimationAdaptive' filesep])
addpath(['DERIVESTsuite' filesep])

%% Load Backup
% if you want to load a backup from a previous estimation (set to [] if
% not)
backup_file = '';

%% Data format
<<<<<<< Updated upstream
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
=======
% Requires a struct array (e.g. data) of size N x 1, where N is number of subjects

% Each entry must contain a structure with at least 3 fields:
% data(n).X : 1xT cell array of T observations.
%    -each entry data(n).X{t} is a matrix of J(t) options x K(t)
%    attributes. This matrix shuld be ordered (along the J dimension) from
%    largest valued alternative to smallest.
%    -assumes that max(J(t)) is the same across all subjects (i.e. each
%    subject sees the biggest choice set)
% data(n).y : Vector of T x 1 Choices
% data(n).J : Vector of T x 1 Choice Set Size

%These additional fields can be populated depending on the model:
%data(n).W : Subject-Specific
%data(n).Z : Attribute
%data(n).K : Number of Attributes per Alternative

%load data/ExampleDataSS
%load data/ExampleData3
%load ~/Dropbox/Projects/Ebbinghaus/code/testdataout.mat
%load ~/Dropbox/Projects/Ebbinghaus/code/subjtestdataout.mat
%load ~/Dropbox/Projects/Ebbinghaus/code/subjtestdataout2.mat
%load ~/Dropbox/Projects/Ebbinghaus/code/realdatabalancedout.mat
%load ~/Dropbox/Projects/Ebbinghaus/code/realdatabalancedout.mat
%load ~/Dropbox/Projects/IIA-NHB/gluth/code/gluthdataout.mat
%d='~/Dropbox/Projects/IIA-NHB/antonio/';
%d='~/Dropbox/Projects/IIA-NHB/kenway/';
d='~/Dropbox/Projects/IIA-NHB/gluth/';

%load([d,'kenwaydataout.mat']);
load([d,'gluthdataout.mat']);
%load([d,'glutheyedataout.mat']);
%load([d,'antoniodataoutHALF.mat']);
%load([d,'antoniodataoutFULL.mat']);
>>>>>>> Stashed changes

%% Estimation Parameters
% You can create your own parameters that will be passed to the likelihood
% function and the particles initialization functions

% Models to use for estimations
% Each model should have a corresponding entry in the following files:

% ProbaChoice : The likelihood of one observation given a model and particle
% setRestrictions : where the parameter space is defined

%{'MNP'}: MNP
%{'DNw'}: allow asymmetric weights
%{'DNwb'}: allow asymmetric weights and beta %is this identified? does not appear to be in some cases....
%{'DN'}: basic DN, omega can be negative
%{'DNb'}: basic DN with beta, omega must be positive
%{'PDNNew'}
%{'RemiStand';'HierarchicalProbit'}
%{'Range'}
<<<<<<< Updated upstream
opts.Models = {'DNb'}

opts.Prob='Probit'; %'Probit','Logit','GHK', 'HP' Is the covariance matrix restricted to be independent?
opts.names={'kappa','sigma','omega','a','b','w2'};
opts.i0=3; %normalize w.r.t. altnerative...
=======
%{'Ebb'}
opts.Models = {'DN'}

opts.Prob='Logit' %'Probit','Logit','GHK', 'HP' Is the covariance matrix restricted to be independent?

opts.names={'kappa','sigma','omega','a','b','wx'};
if ~strcmp(opts.Prob,'Logit')
    opts.names=[opts.names,{'c1','c2'}];
end

>>>>>>> Stashed changes
opts.Hier=logical([0 0 0 0 0 0]); %which parameters to make hierarchical: kappa, sigma, omega, alpha, beta;

%opts.cluster=ones(length(data),1); %all in same cluster (i.e. pooled)
opts.cluster=length(data); % each in own cluster

% opts.Sinz=0;
% opts.weights=0;

%%%%Initial Parameter Values (the number of elements must match the number
%%%%of parameters needed in the model above, +1 for every hierarchical
%%%%parameter)
<<<<<<< Updated upstream
theta0=[];
=======
%theta0=[0.027 0.214 3.458 0.079];
%theta0=[0.012 0.412 513.7 0];
%theta0=[0.012 0.412 25.74];
%theta0=[0 0.5 51 0];
theta0=[0.65754   0  ];

%theta0=[0 0.4923 19 -.032];
%theta0=[0.012 0.412 25.74];
%theta0=[0.012 0.412];
%theta0=[0.0015991      1.7574      1.2482    0.037883  7.7886e-09    0.048639];
>>>>>>> Stashed changes
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
opts.R=200; %number of draws for Halton
opts.numGHKdraws=5000; %number of draws for GHK
opts.ses=0; %calculate ses? Set to zero if you just want estimates (much faster)
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

<<<<<<< Updated upstream
save MLEout

%% Plot Densities from any hierarchical parameters
=======
save([d,'MLEout',opts.Models{:},opts.Prob,'.mat'])
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
>>>>>>> Stashed changes
if any(opts.Hier)
   MLEout.parh(opts.Hier) 
   
   d=gampdf(repmat(0:.01:1,sum(opts.Hier),1), repmat(MLEout.parh(opts.Hier==1)',1,101),repmat(MLEout.parh(end-sum(opts.Hier)+1:end)',1,101));
   %d=logncdf(0:.01:4, MLEout.parh(opts.Hier==1),MLEout.parh(end-sum(opts.Hier)+1:end));
   figure(1)
   plot(0:.01:1,d)
   nanmean(repmat(0:.01:1,sum(opts.Hier),1) .* d)
end
