function [ logLik ] = LogLikelihood( Xs, ChoiceList, subj , model , particleIn, opts, algorithmOpts )
%LIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here

%%THIS FUNCTION IS NOT CURRENTLY USED

T = numel(Xs);
Jt=cellfun(@length,Xs);

data.X = Xs;
data.Mi = {};
for t=1:T
    [data.J(t), data.K(t)] = size(data.X{t});
    temp=eye(Jt(t)-1); 
    for i=1:Jt(t)
        M{i}=[temp(:,1:i-1) -1*ones(Jt(t)-1,1) temp(:,i:Jt(t)-1)];
    end 
    data.Mi{t}=M{ChoiceList(t)}(1:Jt(t)-1,1:Jt(t)); 
end


par(algorithmOpts.LB==algorithmOpts.UB)=algorithmOpts.r0;%Set the restricted variables.
if any(algorithmOpts.LB~=algorithmOpts.UB)
    par(algorithmOpts.LB~=algorithmOpts.UB)=particleIn.theta;%Set the unrestricted variables
end
particleOut.theta = par; %wrap parameters into a structure
Pi=ProbaChoice(data, particleOut, algorithmOpts.model, opts );
logLik = sum(log(Pi));

% %CellfunVersion
% proba_c = @(X) ProbaChoice( X, subj , model , particle, param );
% probas = cellfun(proba_c,Xs,'UniformOutput',0);
% for t = 1:T
%    proba_choice = probas{t};
%    logLik = logLik + log(proba_choice(ChoiceList(t)));
% end


% for t = 1:T
%    proba_choice = ProbaChoice( Xs{t}, subj , model , particle, param );
%    logLik = logLik + log(proba_choice(ChoiceList(t)));
% end

end

