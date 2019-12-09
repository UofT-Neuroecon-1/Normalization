function opts=setRestrictions(modelIn,Jmax,opts)

%Take as argument the model to be estimated. This model can be a
%restricted vesion of a more general model.
%Return 
%modelOut: the function handle to be used to evaluate utilities
%LB: Lower bound on parameters
%UB: Upper bound on parameters

LB.k=1;UB.k=1;
LB.s=0;UB.s=inf;

if strcmp(modelIn,'MNP')
%omega
LB.w=0;UB.w=0;
LB.w2=0;UB.w2=0;
%kappa
LB.a=1;UB.a=1;
LB.b=1;UB.b=1;
opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DN')
%omega
LB.w=0;UB.w=inf;
LB.w2=0;UB.w2=0; 
%kappa
LB.a=1;UB.a=1;
LB.b=1;UB.b=1;
opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNb')
%omega
LB.w=0;UB.w=inf;
LB.w2=0;UB.w2=0;
%kappa
LB.a=1;UB.a=1;
LB.b=0;UB.b=inf;
opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNw')
%omega
LB.w=0;UB.w=inf;
%kappa
LB.a=1;UB.a=1;
LB.b=1;UB.b=1;
LB.w2=0;UB.w2=inf; 
opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'DNwb')
%omega
LB.w=0;UB.w=inf;
%kappa
LB.a=1;UB.a=1;
LB.b=0;UB.b=inf;
LB.w2=-inf;UB.w2=inf; 
opts.modelF='DN'; %use this for estimation
elseif strcmp(modelIn,'Range')
%omega
LB.w=0;UB.w=inf;
LB.w2=1;UB.w2=1;
%kappa
LB.a=1;UB.a=1;
LB.b=1;UB.b=1;
opts.modelF='Range'; %use this for estimation
elseif strcmp(modelIn,'Rangew')
%omega
LB.w=0;UB.w=inf;
LB.w2=-inf;UB.w2=inf;
%kappa
LB.a=1;UB.a=1;
LB.b=1;UB.b=1;
opts.modelF='Range'; %use this for estimation
end


%Define Covariance Matrix
scale=.5; %Set scale of covariance matrix
M=[ -1*ones(Jmax-1,1) eye(Jmax-1)];
L=chol(M*scale*eye(Jmax)*M')';
L=L(L~=0);
cind=L(2:end)';
        
%[1 
% c1 c3
% c2 c4 c5]% 
 

if strcmp(opts.Prob,'Probit') || strcmp(opts.Prob,'Logit') || opts.setsize~=0
    disp('Restricting covariance matrix to be independent.')
    LB.c=cind;
    UB.c=cind;
else
    temp=-inf(Jmax-1);
    temp=tril(temp,-1);
    temp=temp(L~=0)';
    
    %LB.c=[-inf, 0];
    LB.c=temp(2:end);
    UB.c=inf(1,length(cind));
    
    myseed = 20110710; %Set seed for GHK
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',myseed));
    E=rand(Jmax-1,opts.R); %Draw    
end

if any(opts.Hier) %set scale hyperparameter(s) to be between 0 and inf. The mean or shape will be passed in the respective loation LB.par and UB.par
    LB.Gamma=zeros(1,sum(opts.Hier));
    UB.Gamma=inf(1,sum(opts.Hier));
else
    LB.Gamma=[];
    UB.Gamma=[];
end


opts.LB=[LB.k LB.s LB.w LB.a LB.b LB.w2 LB.c LB.Gamma];
opts.UB=[UB.k UB.s UB.w UB.a UB.b UB.w2 UB.c UB.Gamma];
     
end