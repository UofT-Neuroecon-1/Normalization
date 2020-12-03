function opts=setRestrictions(modelIn,Jmax,opts)

%Take as argument the model to be estimated. This model can be a
%restricted vesion of a more general model.
%Return 
%modelOut: the function handle to be used to evaluate utilities
%LB: Lower bound on parameters
%UB: Upper bound on parameters



if strcmp(modelIn,'MNP')
    opts.toEst={'sigma'};
    LB=0;
    UB=Inf;
    opts.toRestr={'omega','alpha','beta','wx'};
    thetaR=[0 1 1 0];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'Logit')
    opts.toEst={'sigma'};
    LB=0;
    UB=Inf;
    opts.toRestr={'omega','alpha','beta','wx'};
    thetaR=[0 1 1 0];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DN')
    opts.toEst={'sigma','omega'};
    LB=[0 -Inf];
    UB=[Inf Inf];
    opts.toRestr={'alpha','beta','wx'};
    thetaR=[1 1 0];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNw3')
    opts.toEst={'sigma','omega1','omega2','omega3'};
    LB=[0 -Inf -Inf -Inf];
    UB=[Inf Inf Inf Inf];
    opts.toRestr={'alpha','beta'};
    thetaR=[1 1];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNw2')
    opts.toEst={'sigma','omega1','omega2'};
    LB=[0 -Inf -Inf ];
    UB=[Inf Inf Inf ];
    opts.toRestr={'alpha','beta','omega3'};
    thetaR=[1 1 0];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNb')
    opts.toEst={'sigma','omega','beta'};
    LB=[0 0 0];
    UB=[Inf Inf 5000];
    opts.toRestr={'alpha','wx'};
    thetaR=[1 0];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNw')
    opts.toEst={'sigma','omega','wx'};
    LB=[0 -Inf 0];
    UB=[Inf Inf Inf];
    opts.toRestr={'alpha','beta'};
    thetaR=[1 1];
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNwb')
    opts.toEst={'sigma','omega','beta','wx'};
    LB=[0 -Inf 0 0];
    UB=[Inf Inf 5000 Inf];
    opts.toRestr={'alpha'};
    thetaR=1;
    opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'Range')
    opts.toEst={'sigma','omega'};
    LB=[0 -Inf];
    UB=[Inf Inf];
    opts.toRestr={'alpha','beta','wx'};
    thetaR=[1 1 1];
    opts.modelF='Range'; %use this for estimation
elseif strcmp(modelIn,'Rangew')
    opts.toEst={'sigma','omega','wx'};
    LB=[0 -Inf -Inf];
    UB=[Inf Inf Inf];
    opts.toRestr={'alpha','beta'};
    thetaR=[1 1];
    opts.modelF='Range'; %use this for estimation
elseif strcmp(modelIn,'Ebb')
    opts.toEst={'sigma','omega','wx'};
    LB=[0 -Inf 0];
    UB=[Inf Inf Inf];
    opts.toRestr={'alpha','beta','webb'};
    thetaR=[1 1 0];
    opts.modelF='Ebb'; %use this for estimation
elseif strcmp(modelIn,'Ebb2')
    opts.toEst={'sigma','omega','wx','webb'};
    LB=[0 -Inf 0 -Inf];
    UB=[Inf Inf Inf Inf];
    opts.toRestr={'alpha','beta'};
    thetaR=[1 1];
    opts.modelF='Ebb'; %use this for estimation
end
 
if Jmax==2
    opts.toRestr=[opts.toRestr {'c1'}];
else
    %Define Covariance Matrix
    M=[ -1*ones(Jmax-1,1) eye(Jmax-1)];
    L=chol(M*opts.scale*eye(Jmax)*M')';
    L=L(L~=0);
    cind=L(2:end)';

    %[1 
    % c1 c3
    % c2 c4 c5]% 
    if strcmp(opts.Prob,'Probit') || opts.setsize~=0
        disp('Restricting covariance matrix to be independent.')
        thetaR=[thetaR cind];
        
        for j=1:length(cind)
            opts.toRestr=[opts.toRestr {['c' num2str(j)]}];
        end
        
    elseif strcmp(opts.Prob,'GHK')
        temp=-inf(Jmax-1);
        temp=tril(temp,-1);
        temp=temp(L~=0)';


        for j=1:Jmax-1
            opts.toEst=[opts.toEst {['c' num2str(j)]}];
        end

        LB=[LB temp(2:end)];
        UB=[UB inf(1,length(cind))];
%         myseed = 20110710; %Set seed for GHK
%         RandStream.setGlobalStream(RandStream('mt19937ar','seed',myseed));
%         opts.GHKdraws=rand(Jmax-1,opts.numGHKdraws); %Draw
        
        rng(1,'twister') % set random number seed
        p = haltonset((Jmax-1),'Skip',1e3,'Leap',1e2);  %Halton Sequence
        opts.GHKdraws=net(p,opts.numGHKdraws)';  
    end
end

if ~isempty(opts.Hier) %set scale hyperparameter(s) to be between 0 and inf. The mean or shape will be passed in its respective loation LB.par and UB.par
    LB=[LB zeros(1,length(opts.Hier))];
    UB=[UB inf(1,length(opts.Hier))];
    opts.scale=opts.scale; %set scale lower so that density on omega can evalute away from zero if need be
    %set opts.scale=opts.scale/50;
    opts.Prob='Logit';
    disp('Hierarchical Model: Setting error distribution to Logit')%code for other ditributions might not be vectorized.
end

opts.LB=LB;
opts.UB=UB;
opts.thetaR=thetaR;
     
end