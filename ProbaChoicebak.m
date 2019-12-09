function [ proba_choice ] = ProbaChoice( X, subj , model , particle, param )
% X: J x K matrix of choice set
% model: the true model to use
% (returns) proba_choice : Jx1 vector of choice probabilities
% it is recommended to make sure that the probability of a choice is not 0

J = size(X,1);
K = size(X,2);

if strcmp(model,'Logit')
    %True params
    alpha = particle.theta(subj,1);
    Beta = (param.attrSign .* particle.theta(subj,2:end))';
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    for j=1:J
        v(j) = sum(Beta .* u_x(j,:)'); 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
    %mixture 99.9% model and 0.1% unif
    proba_choice = 0.99 .* proba_choice + 0.01/J;
elseif strcmp(model,'PDNNew')
    %True params
    alpha = particle.theta(subj,1);
    sigma = particle.theta(subj,2);
    omega = particle.theta(subj,3:3+K-1);
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u = (param.attrSign .* u_x)';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y)),1);%./(J-1);
        v(j) = norm_coefs * unnorm_u(:,j); 
    end
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
    %mixture 99.9% model and 0.1% unif
    proba_choice = 0.99 .* proba_choice + 0.01/J;
 elseif strcmp(model,'PDNProbit')
    %True params
    alpha = particle.theta(subj,1);
    sigma = particle.theta(subj,2);
    omega = particle.theta(subj,3:3+K-1);
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u = (param.attrSign .* u_x)';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y)),1);%./(J-1);
        v(j) = norm_coefs * unnorm_u(:,j); 
    end
    
        temp=eye(J-1); 
    for i=1:J
        M{i}=[temp(:,1:i-1) -1*ones(J-1,1) temp(:,i:J-1)];
    end

     Mi=M{y}(1:J-1,1:J); 

    proba_choice(y) = calcPiInd(Mi,v,J);
    %mixture 99.9% model and 0.1% unif
    proba_choice(y) = 0.99 .* proba_choice(y) + 0.01/J;
elseif strcmp(model,'DNv')
    %True params
    J = size(X{1},1);
    K = size(X{1},2);

    alpha = particle.theta(subj,1);
    sigma = particle.theta(subj,2);
    omega = particle.theta(subj,3:3+K-1);

    f = @(x) (x.^alpha);
    denom=@(x) (sigma + omega*sum(x) );
    sumv=cellfun(denom,X,'uniformoutput',false);
    v=cellfun(@rdivide,cellfun(f,X,'uniformoutput',false),sumv,'uniformoutput',false);

    proba_choice = calcPiInd(y,v,J); %y is Mi
    %mixture 99.9% model and 0.1% unif
    proba_choice = 0.99 .* proba_choice + 0.01/J;

elseif strcmp(model,'range')
        
        denom=@(x) (sigma + omega*(max(x)-min(x)));
    
elseif strcmp(model,'RemiStand')
    %True params
    alpha = particle.theta(subj,1);
    sigma =particle.theta(subj,2);
    Beta = (param.attrSign .* particle.theta(subj,3:3+K-1))';
    %utility computation
    x_mean = mean(X,1);
    sd_x = std(X,[],1) * alpha + sigma;
    x_standardized = (X - x_mean) ./ sd_x;
    attr_signal = 1 ./ (1+exp(-x_standardized));
    v = attr_signal * Beta;
    v = v - max(v); %avoid overflow
    sum_exp_v = sum(exp(v));
    proba_choice = exp(v)./sum_exp_v;
    %mixture 99.9% model and 0.1% unif
    proba_choice = 0.99 .* proba_choice + 0.01/J;
elseif strcmp(model,'HierarchicalProbit')
    % [alpha sigma Omega(1,K)]
    %True params
    alpha = particle.theta(subj,1);
    sigma = particle.theta(subj,2);
    Beta = (param.attrSign)';
    omega = particle.theta(subj,3:3+K-1);
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u = (param.attrSign .* u_x)';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y) ),1);%./(J-1);
        v(j) = norm_coefs * unnorm_u(:,j); 
    end
    v = v - max(v); %avoid overflow
    VarCov = eye(J,J);
    %Cholesky decomp + check positive def
    [CholeskyUpper,psd] = chol(VarCov);
    while psd
        VarCov = VarCov + 0.00001 * eye(size(VarCov,1));
        [CholeskyUpper,psd] = chol(VarCov);
    end
    sim = repmat(v',[size(param.NormDraw,1) 1]) + param.NormDraw(:,1:J) * CholeskyUpper;
    proba_choice = mean(sim == max(sim,[],2));
    %mixture 99.9% model and 0.1% unif
    proba_choice = 0.99 .* proba_choice + 0.01/J;
else
    error('ProbaChoice : unknown model');
end

function Pi=calcPiInd(Mi,v,J)

       T=size(v,2);

       
       if T==1
            vi=Mi*v;
            [x, w]=GaussHermite(100);

         %   vi = cell2mat(viC);
            zz2=bsxfun(@minus,-sqrt(2).*vi,repmat(-sqrt(2)*reshape(x,[1 100]), [J-1,1]));
            aa2=prod(normcdf(zz2),1);
            Pi=sum(bsxfun(@times,w',squeeze(aa2)),2)./sqrt(pi);
       else
        
            viC = cellfun(@mtimes, Mi, v, 'UniformOutput', false); %cell version

            [x, w]=GaussHermite(100);

            vi = cell2mat(viC);
            zz2=bsxfun(@minus,-sqrt(2).*vi,repmat(-sqrt(2)*reshape(x,[1 1 100]), [J-1,T,1]));
            aa2=prod(normcdf(zz2),1);
            Pi=sum(bsxfun(@times,w',squeeze(aa2)),2)./sqrt(pi);
       end
end


end

