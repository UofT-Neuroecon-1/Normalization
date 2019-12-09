function [ Pi ] = ProbaChoice( data, particle,opts )
    % data.X: J x K matrix of choice set
    % data.y: Tx1 of choices
    % model: the model to use
    % (returns) Pi : Jx1 vector of choice probabilities
    %F: function handle for model
    

    J = data.J;
    K = data.K;
    
    F=eval(['@' opts.modelF]);
    P=eval(['@calcPi' opts.Prob]); %Construct function handle for probability function
        
%       if opts.indep %%%%Gaussian Quadrature       
%             %If # of alts not constant over sample, calc Pi's individually
%             if opts.SS==1 %set size changes on each trial, so run for loop over trials (might be able to speed up like below...)
%                 Pi=calcPiSS(Mi,Z);
%             else %set size is constant
%                 Pi=calcPiInd(Mi,Z,J);
%             end                  
%        else %Use GHK simulation
%             if opts.SS==1
%                 Pi=calcPiGHKC(Mi,Z,J,c,E);
%             else
%                 Pi=calcPiGHK(y,cell2mat(Z)',c,E);
%             end
%        end
% 
    
    Pi=F(data.X,particle.theta); %Get Probs
 
    Pi(Pi==0)=realmin;  % make sure that the probability of a choice is not 0
    

    function Pi=Logit(X,theta)
        %True params
        alpha = theta(1);
        Beta = (opts.attrSign .* theta(2:end))';
        %utility computation
        u_x = X.^alpha;
        v = zeros(J,1);
        for j=1:J
            v(j) = sum(Beta .* u_x(j,:)'); 
        end
        v = v - max(v); %avoid overflow
        sum_exp_v = sum(exp(v));
        Pi = exp(v)./sum_exp_v;
        %mixture 99.9% model and 0.1% unif
        Pi = 0.99 .* Pi + 0.01/J;
    end

    function Pi=PDNNew(X,theta)
        %True params
        alpha = theta(1);
        sigma = theta(2);
        omega = theta(3:3+K-1);
        %utility computation
        u_x = X.^alpha;
        v = zeros(J,1);
        unnorm_u = (opts.attrSign .* u_x)';
        for j=1:J
            u_y = u_x;
            u_y(j,:)=[];
            norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y)),1);%./(J-1);
            v(j) = norm_coefs * unnorm_u(:,j); 
        end
        v = v - max(v); %avoid overflow
        sum_exp_v = sum(exp(v));
        Pi = exp(v)./sum_exp_v;
        %mixture 99.9% model and 0.1% unif
        %proba_choice = 0.99 .* proba_choice + 0.01/J;
    end

    function Pi=PDNProbit(X,theta)
        %True params
        alpha = theta(1);
        sigma = theta(2);
        omega = theta(3:3+K-1);
        %utility computation
        u_x = X.^alpha;
        v = zeros(J,1);
        unnorm_u = (opts.attrSign .* u_x)';
        for j=1:J
            u_y = u_x;
            u_y(j,:)=[];
            norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y)),1);%./(J-1);
            v(j) = norm_coefs * unnorm_u(:,j); 
        end

        Pi(data.y) = P(data.Mi,v,J);
        %mixture 99.9% model and 0.1% unif
        %proba_choice(y) = 0.99 .* proba_choice(y) + 0.01/J;
    end

    function Pi=DN(X,theta)
        
        par=opts.LB;
        if any(opts.LB(1:8)~=opts.UB(1:8))
            par(opts.LB(1:8)~=opts.UB(1:8))=theta;%Set the unrestricted variables to be those passed to the function.
        end

        s = par(2);
        w = par(3); % all weights
        a = par(4);
        b = par(5);
        w2= par(6); %additional own weight
        
        if b<0
            b=0.0000000001; %hack so that s.e. calculation (hessian()) can send negative b. In estimation, beta is restricted to be >0 so doesn't matter.
        end

        f = @(x) (x.^a);
        %denom=@(x) (s + w*sum(x) );
        
        %denom=@(x) (s + (w + w2*eye(length(x)))* x );

        denom=@(x) (s + vecnorm((w + w2*eye(length(x))).*x',b,2) ); %Note that w2 is the "extra" weight relative w. So in reporting results, add them together for w_i
        %denom=@(x) (s + ((w + w2*eye(length(x))).*x.^b').^(1/b) );
        
        %denom=@(x) (sigma + omega*norm(x,b) );
        %vecnorm(cell2mat(X)',2)
        sumv=cellfun(denom,X,'uniformoutput',false);
        v=cellfun(@rdivide,cellfun(f,X,'uniformoutput',false),sumv,'uniformoutput',false);

        Pi = P(data.Mi,v,data.J); %y is Mi
    end

    function Pi=Range(X,theta)
 
        par=opts.LB;
        if any(opts.LB(1:8)~=opts.UB(1:8))
            par(opts.LB(1:8)~=opts.UB(1:8))=theta;%Set the unrestricted variables to be those passed to the function.
        end

        s = par(2);
        w = par(3); % all weights
        a = par(4);
        b = par(5);
        w2= par(6); %additional own weight
        
%             kappa=par(1:Q);
%             sigma=par(Q+1);
%             omega=par(Q+2);                
%             a=par(Q+3);
%             b=par(Q+4);
%             w2=par(Q+5);
%             c=par(Q+6:end);
% 
%             %1  0   ... 0   
%             %c1 c2  ... 0
%             %c3 ... ... cend

        f = @(x) (x.^a);
        %denom=@(x) (s + w*sum(x) );
        
        denom=@(x) (s + w*(max(x)-w2*min(x)));
        
        %denom=@(x) (sigma + omega*norm(x,b) );
        %vecnorm(cell2mat(X)',2)
        sumv=cellfun(denom,X,'uniformoutput',false);
        v=cellfun(@rdivide,cellfun(f,X,'uniformoutput',false),sumv,'uniformoutput',false);

        Pi = P(data.Mi,v,data.J); %y is Mi
            
    end

    function Pi=RemiStand(X,theta)
        %True params
        alpha =theta(1);
        sigma =theta(2);
        Beta = (opts.attrSign .* theta(subj,3:3+K-1))';
        %utility computation
        x_mean = mean(X,1);
        sd_x = std(X,[],1) * alpha + sigma;
        x_standardized = (X - x_mean) ./ sd_x;
        attr_signal = 1 ./ (1+exp(-x_standardized));
        v = attr_signal * Beta;
        v = v - max(v); %avoid overflow
        sum_exp_v = sum(exp(v));
        Pi = exp(v)./sum_exp_v;
        %mixture 99.9% model and 0.1% unif
        Pi = 0.99 .* Pi + 0.01/J;
    end

    function Pi=HierarchicalProbit(X,theta)
        % [alpha sigma Omega(1,K)]
        %True params
        alpha = theta(1);
        sigma = theta(2);
        Beta = (opts.attrSign)';
        omega = theta(subj,3:3+K-1);
        %utility computation
        u_x = X.^alpha;
        v = zeros(J,1);
        unnorm_u = (opts.attrSign .* u_x)';
        for j=1:J
            u_y = u_x;
            u_y(j,:)=[];
            norm_coefs = sum(1 ./ (sigma + omega .* (u_x(j,:) + u_y) ),1);%./(J-1);
            v(j) = norm_coefs * unnorm_u(:,j); 
        end
       
        Pi = P(v,J,opts.NormDraw);
    end
end

function Pi=calcPiProbit(Mi,v,J)

    T=size(v,2);
    [x, w]=GaussHermite(100);
       
    if T==1 %one trial at a time
        if isa(v,'cell')   
           v = cell2mat(v);
           Mi=cell2mat(Mi);
        end
        vi=Mi*v;

        %zz=bsxfun(@minus,-sqrt(2).*vi,repmat(-sqrt(2)*reshape(x,[1 100]), J-1,1));
        zz=(-sqrt(2).*vi) - reshape(x,[1,100]);
        aa=prod(normcdf(zz),1);
        Pi=sum(bsxfun(@times,w',squeeze(aa)),2)./sqrt(pi);
        
    elseif all(J==J(1)) % all trials together (much faster)   
        viC = cellfun(@mtimes, Mi, v, 'UniformOutput', false); %cell version
        vi = cell2mat(viC);

        %zz=bsxfun(@minus,-sqrt(2).*vi,repmat(-sqrt(2)*reshape(x,[1 1 100]), [J-1,T,1]));
        zz=(-sqrt(2).*vi) - reshape(-sqrt(2)*x,[1,1,100]);
        aa=prod(normcdf(zz));
        Pi=sum(bsxfun(@times,w',squeeze(aa)),2)./sqrt(pi);

    else %all trials, but works with cells in case J is different over trials
        viC = cellfun(@mtimes, Mi, v, 'UniformOutput', false);

        f=@(y) (bsxfun(@minus, -sqrt(2).*y , -sqrt(2)*x'));
        zz=cellfun(f,viC,'uniformoutput',false);
        aa=cellfun(@(x) prod(x,1),cellfun(@normcdf,zz,'uniformoutput',false),'uniformoutput',false);
        Pi=sum(bsxfun(@times,w',cell2mat(aa')),2)./sqrt(pi);

    end  
        %mixture 99.9% model and 0.1% unif
        %Pi = 0.99 .* Pi + 0.01/J;
end

function Pi=calcPiLogit(Mi,v,~)
        vi=cellfun(@mtimes, Mi, v, 'UniformOutput', false);
        sum_exp_v = sum(exp(cell2mat(vi)));
        Pi = 1./(1+sum_exp_v);
end
            

function Pi=calcPiHP(v,J,eps)
        v = v - max(v); %avoid overflow
        VarCov = eye(J,J);
        %Cholesky decomp + check positive def
        [CholeskyUpper,psd] = chol(VarCov);
        while psd
            VarCov = VarCov + 0.00001 * eye(size(VarCov,1));
            [CholeskyUpper,psd] = chol(VarCov);
        end
        sim = repmat(v',[size(eps,1) 1]) + eps(:,1:J) * CholeskyUpper;
        Pi = mean(sim == max(sim,[],2));
        %mixture 99.9% model and 0.1% unif
        Pi = 0.99 .* Pi + 0.01/J;
end

% 
% function Pi=calcPiSS(Mi,Z)  %do it at the cell level because J changes each trial
% 
%         ZiC = cellfun(@mtimes, Mi, Z, 'UniformOutput', false);
%         
%         [x, w]=GaussHermite(100);
%         
%         f=@(y) (bsxfun(@minus, -sqrt(2).*y , -sqrt(2)*x'));
%         zz=cellfun(f,ZiC,'uniformoutput',false);
%         aa=cellfun(@(x) prod(x,1),cellfun(@normcdf,zz,'uniformoutput',false),'uniformoutput',false);
%         Pi=sum(bsxfun(@times,w',cell2mat(aa')),2)./sqrt(pi);
% end
% 
% % function Pi=calcPiSS2(Mi,Z)  %do it at the cell level because J changes each trial
% % 
% %         ZiC = cellfun(@mtimes, Mi, Z, 'UniformOutput', false);
% %         
% %         
% %         Pi=cellfun(@mvncdf,ZiC,'UniformOutput',true)'; This is wrong
% % 
% % end
% 
function Pi=calcPiGHKC(Mi,Z,J,c,E)
        
        %[T,J]=size(Z);
        
%         L = tril(ones(J),0);
%         L(:,1)=zeros(J,1);
%         L(~~L)=[1 c];
%         Omega=L*L';
        
        L1 = tril(ones(J-1),0);
        L1(~~L1)=[1 c]; 
        L=[zeros(J-1,1) L1];
        L=[zeros(1,J); L];
        Omega=L*L';
        
        Vi=cellfun(@mtimes,Mi,Z,'uniformoutput',false); %cellify
        
        %Reset seed inside each evaluation of LL so that no noise inGHK over evaluations 

        cholfun=@(X) (chol(X*Omega*X')');

        Li=cellfun(cholfun,Mi,'uniformoutput',false);

        ghkfun=@(V,L) (ghktrain(V,L,E));
        Pi=cellfun(ghkfun,Vi,Li);  %cellify
end

% function Pi=calcPiGHK(y,Z,c,E)
%         
%         [T,J]=size(Z);
% 
%         L1 = tril(ones(J-1),0);
%         L1(~~L1)=[1 c]; 
%         L=[zeros(J-1,1) L1];
%         L=[zeros(1,J); L];
%         Omega=L*L';
%         
%         
%         %Reset seed inside each evaluation of LL so that no noise inGHK over evaluations 
%         myseed = 20110710;
% 
%         RandStream.setGlobalStream(RandStream('mt19937ar','seed',myseed));
% 
%         for t=1:T %do not parfor the GHK unless you check that the RNG is the same for each t on each call
%             j=y(t);
% 
%             temp=eye(J-1); 
%             Mi=[temp(:,1:j-1) -1*ones(J-1,1) temp(:,j:J-1)];
% 
%             Vi=Mi*Z(t,:)';
%             %Vi=zeros(JJ-1,1);
% 
%                Li=chol(Mi*Omega(1:J,1:J)*Mi')';
% 
% 
%             Pi(t)=ghktrain(Vi,Li,E);   %Sim Pi 
% 
% %                         temp=eye(J-1); 
% %                         Mi=[temp(:,1:i-1) -1*ones(J-1,1) temp(:,i:J-1)];
% 
%             %Pi(t)=GHK(Li,Vi',R);
%         end
% end