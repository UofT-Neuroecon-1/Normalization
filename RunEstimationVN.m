%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation / Model-comparison Example File %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

diary on
addpath(['..' filesep 'EstimationAdaptive' filesep])
addpath(['DERIVESTsuite' filesep])

%% Load Backup
% if you want to load a backup from a previous estimation (set to [] if
% not)
backup_file = '';

%% Data format
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
load data/ExpDataSS
%load data/ExpData3

if ~isfield(data,'Z')
    data(1).('Z') = [];
end
if ~isfield(data,'W')
    data(1).('W') = [];
end
%% Estimation Parameters
% You can create your own parameters that will be passed to the likelihood
% function and the particles initialization functions

% Models to use for estimations
% Each model should have a corresponding entry in the following files:

% ProbaChoice : The likelihood of one observation given a model and particle
% setRestrictions : where the parameter space is defined

%{'MNP'}: MNP
%{'DNw'}: allow asymmetric weights
%{'DNwb'}: allow asymmetric weights and beta %is this identified? does not appear to be in some cases.... (i.e. of beta large, or w<0)
%{'DN'}: basic DN, omega can be negative
%{'DNb'}: basic DN with beta, omega must be positive
%{'PDNNew'}
%{'RemiStand';'HierarchicalProbit'}
%{'Range'}
%{'Ebb'}
%{'Logit'}
%'DNw3'
opts.Models = {'DN'};

opts.Prob='Probit'; %'Probit','Logit','GHK', 'HP' Is the covariance matrix restricted to be independent?

opts.WithinSubject=0;
opts.Hier={}; %which parameters to make hierarchical: kappa, sigma, omega, alpha, beta;
%opts.Hier={'sigma','omega'};
opts.HierDist={'gamma','normal'};

opts.cluster=length(data); % each in own cluster

%%%%Initial Parameter Values (the number of elements must match the number
%%%%of parameters needed in the model above, +1 for every hierarchical
%%%%parameter)

theta0=[.783    0.027];

%%%%Estimation Specific Parameters
%for particle filter
opts.G = 3; % Number of particles group
opts.P = 128; % Number of particles per group
opts.Adaptive = true; % Use the adaptive SMC (see Durham, Geweke 2014).
opts.Msteps = 10; % Number of mutate steps
opts.ress_threshold = 0.8;

%for ML
opts.getP = false; %Set this to true if you just want to evaluate at initial parameter vector
opts.numInit = 1; %1: run using gradient- method first at theta0. If >1, use random starting points with gradient free method.
opts.R = 50000; %number of draws for hierchical
opts.numGHKdraws = 5000; %number of draws for GHK
opts.ses = true; %calculate ses? Set to zero if you just want estimates (much faster)
opts.Display='iter-detailed';

%%%% Model specific parameters
%opts.NormDraw = mvnrnd(zeros(4,1),eye(4),1000); % Pre-draw errors for GHK (if needed)

opts.attrVals{1}= (0:4);          opts.attrNames{1}="Bid";
opts.attrSign = [1];
opts.scale=0.5; %set scale/std.dev. of error term to 1/2 so that differenced is normalized to 1. This is done in sim code as well
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

%% Max Likelihood 
% Maximum Likelihood (Within)
if opts.WithinSubject
    for s = 1:numel(data)
        MLEout(s) = MLestimation(data(s),theta0,opts);

        save([d,'MLEout',opts.Models{:},opts.Prob,'WithinSubject','.mat'])
    end
else
    %Maximum Likelihood (Pooled or Hierarchical)
    MLEout = MLestimation(data,theta0,opts);

    save([d,'MLEout',opts.Models{:},opts.Prob,'.mat'])
end
%%

if opts.WithinSubject
    for s=1:numel(data)
        sigmas(s)=MLEout(s).parh(1);
        omegas(s)=MLEout(s).parh(2);
    end

figure(1)
clf
subplot(2,1,1)
bar(sigmas)
ylabel('sigma')
subplot(2,1,2)
bar(omegas)
ylabel('omega')
sum(omegas>0)
mean(omegas)
end
%%
if ~isempty(opts.Hier)
   MLEout.parh 
   
for k=1:length(opts.Hier)  
    parh = MLEout.parh(strcmp(opts.Hier(k), [MLEout.toEst,opts.Hier]));
    if strcmp(opts.HierDist(k),'gamma')
            density(k,:)=gampdf(0:1:100, parh(1),parh(2));
        if strcmp(opts.Hier(k),'omega')
           DNvariance=gampdf([0:1:100]./8.2036, parh(1),parh(2));
        end
        parmean(k)=2*mean([1:1:100] .* gampdf([1:1:100],parh(1),parh(2))) %scale by 2 to match default scale in model
    elseif strcmp(opts.HierDist(k),'normal')
            density(k,:)=normpdf(-5:.01:5, parh(1),parh(2));
        if strcmp(opts.Hier(k),'omega')
           DNvariance=gampdf([-5:.01:5]./8.2036, parh(1),parh(2));
        end
        parmean(k)=2*mean([-5:.01:5] .* normpdf([-5:.01:5],parh(1),parh(2))) %scale by 2 to match default scale in model
    end
end
   figure(1)
   clf
   plot([-5:.01:5],density,'linewidth',2) 
   %set(gca,'XTicks',0:1:100);
   xticks(-5:1:5)
   ylim([0,10])
   set(gca,'XTickLabels',[-5:1:5]*opts.scale*2); %x2 because default scale for all parameters is 0.5
    hold on
    for k=1:length(opts.Hier)
        plot([parmean(k)/opts.scale/2 parmean(k)/opts.scale/2], ylim, 'k--') %/2 because default scale for all parameters is 0.5
    end
    text(0,7,['Mean: ',num2str(parh(1))])
    text(0,6,['Variance: ',num2str(parh(2))])
   hold off
   box off
   ylabel('Probability Density','interpreter','latex')
   
   legend(cellfun(@(x) ['$\' x '$'],opts.Hier,'UniformOutput',false),'interpreter','latex')
   
   figure(2)
   clf
   k=find(strcmp(opts.Hier,'omega'));
   variance=density;
   variance(k,:)=DNvariance;
   
   if strcmp(opts.HierDist(k),'gamma')
       plot([0:1:100],variance,'linewidth',2) 
       %set(gca,'XTicks',0:1:100);
       xticks(0:10:100)
       set(gca,'XTickLabels',[0:10:100]*opts.scale*2); %x2 because default scale for all parameters is 0.5
       ylim([0,.1])
   elseif strcmp(opts.HierDist(k),'normal')
       plot([-5:.01:5],variance,'linewidth',2) 
       %set(gca,'XTicks',0:1:100);
       xticks(-5:1:5)
       set(gca,'XTickLabels',[-5:1:5]*opts.scale*2); %x2 because default scale for all parameters is 0.5
   end
   hold on
   set(gca,'ColorOrderIndex',1)
    plot([parmean(1)/opts.scale/2 parmean(1)/opts.scale/2], ylim, '--') %/2 because default scale for all parameters is 0.5
    plot([parmean(2)/opts.scale/2*8.2036 parmean(2)/opts.scale/2*8.2036], ylim, '--') %/2 because default scale for all parameters is 0.5
   hold off
   box off
   ylabel('Probability Density of Variance Contribution','interpreter','latex')
   
   legend(cellfun(@(x) ['$\' x '$'],opts.Hier,'UniformOutput',false),'interpreter','latex')
end

diary off
