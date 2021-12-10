function [ logPi ] = ProbaChoice( data, particle,opts )
    % data.X: J x K matrix of choice set
    % data.y: Tx1 of choices
    % model: the model to use
    % (returns) Pi : Jx1 vector of choice probabilities
    %opts.modelF: contains function handle for model
    

    J = data.J;
    K = data.K;
    

    
    par=[particle.theta repmat(opts.thetaR,1,1,size(particle.theta,3))];
    names=[opts.toEst opts.toRestr];
    
    %Which model to use?
    NormalizationFunction=eval(['@' opts.modelF]);
    
    %Get Valuations
    v=NormalizationFunction(data,par,names,opts);
   
    %Calculate Choice Probabilities
    switch opts.Prob
        case 'GHK'
            if opts.setsize==1
                Pi=calcPiGHKC(data.Mi,v,data.J,par(find(startsWith(names, 'c'))));
                logPi=log(Pi);
            else
                Pi=calcPiGHK(data.y,v',par(find(startsWith(names, 'c'))),opts.GHKdraws);
                logPi=log(Pi);
            end
        case 'Linear'
            logPi=calcPiLinear(data.y,v');
        otherwise               
            P=eval(['@calcPi' opts.Prob]); %Construct function handle for probability function based on label in opts.Prob
            Pi = P(data.Mi,v,data.J,opts.scale); %y is Mi
            logPi=log(Pi);
    end
        

    function Pi=Range(data,theta)
 
        par=[theta opts.thetaR];
        names=[opts.toEst opts.toRestr];
        
        s = par(find(strcmp(names, 'sigma')));
        w = par(find(strcmp(names, 'omega'))); % all weights
        a = par(find(strcmp(names, 'alpha')));
        b = par(find(strcmp(names, 'beta')));
        wx= par(find(strcmp(names, 'wx'))); %additional own weight
        
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
        
        denom=@(x) (s + w*(max(x)-wx*min(x)));
        
        %denom=@(x) (sigma + omega*norm(x,b) );
        %vecnorm(cell2mat(X)',2)
        sumv=cellfun(denom,data.X,'uniformoutput',false);
        numerator=cellfun(f,data.X,'uniformoutput',false);
        v=cellfun(@rdivide,numerator,sumv,'uniformoutput',false);
        
%         P=eval(['@calcPi' opts.Prob]);
%         Pi = P(data.Mi,v,data.J); %y is Mi
%         logPi=log(Pi);
            
    end

function Pi=Ebb(X,prob,ebbc,theta)
        
        %theta is a vector, it could be made a matrix to speed up
        %calculations of hierarchical effects
    
        par=[theta opts.thetaR];
        names=[opts.toEst opts.toRestr];
        
        s = par(find(strcmp(names, 'sigma')));
        w = par(find(strcmp(names, 'omega'))); % all weights
        a = par(find(strcmp(names, 'alpha')));
        b = par(find(strcmp(names, 'beta')));
        wx= par(find(strcmp(names, 'wx'))); %additional own weight
        webb= par(find(strcmp(names, 'webb'))); %additional own weight
        
       
        if b<0
            b=0.0000000001; %hack so that s.e. calculation (hessian()) can send negative b. In estimation, beta is restricted to be >0 so doesn't matter.
        end

        f = @(x) (x.^a);
        %denom=@(x) (s + w*sum(x) );
        
        denom=@(x,ebb) ((s) + ((w+ webb*ebb) + wx*eye(length(x)))* x );
        %denom=@(x,ebb) (s + ((w+ webb*ebb) + w2*eye(length(x)))* x );
        %denom=@(x,ebb) (s + (w + w2*eye(length(x)))* x );
        %denom=@(x,ebb) (s + (w2*eye(length(x)))*w*x );
        
        %denom=@(x) (s + vecnorm((w + w2*eye(length(x))).*x',b,2) ); %Note that w2 is the "extra" weight relative w. So in reporting results, add them together for w_i
        %denom=@(x) (s + ((w + w2*eye(length(x))).*x.^b').^(1/b) );
        
        %denom=@(x) (sigma + omega*norm(x,b) );

        sumv=cellfun(denom,X,ebbc,'uniformoutput',false);
        numerator=cellfun(f,X,'uniformoutput',false);
        v=cellfun(@rdivide,numerator,sumv,'uniformoutput',false);
        eu=cellfun(@times,prob,v,'uniformoutput',false);
        
%         P=eval(['@calcPi' opts.Prob]);
%         Pi = P(data.Mi,eu,data.J); %y is Mi
%         logPi=log(Pi);
    end
end

%Get Choice Probs given the model above

function Pi=calcPiProbit(Mi,v,J,scale)

    T=size(v,2);
    
    
    if size(v,1)==2 %Binary Choice
        if all(J==J(1)) % all choice sets same size, then all trials together (much faster)   
            vi=pagemtimes(Mi, permute(v,[1,3,2]));
            %vi=cellfun(@mtimes, Mi, v)';
            Pi=normcdf(-vi./scale);
        else
            error('Different size choice sets. This Probit code doesnt work for binary choice, use Logit')
        end
        
    else %Choice sets larger than binary, use quadrature
    
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
%             viC = cellfun(@mtimes, Mi, v, 'UniformOutput', false); %cell version
%             vi = cell2mat(viC);
            vi=squeeze(pagemtimes(Mi, permute(v,[1,3,2])));
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
end

function Pi=calcPiLogit(Mi,v,~,scale)
        %vi=cellfun(@mtimes, Mi, v, 'UniformOutput', false);
        %cellfun(@exp,vi','UniformOutput', false)
        
        vi=pagemtimes(Mi, permute(v,[1,3,2]));
        sum_exp_v = sum(exp(vi./scale),1);
        %sum_exp_v = cellfun(@(x) sum(exp(x./scale),1),vi,'UniformOutput', false);
        Pi = 1./(1+squeeze(sum_exp_v)');
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

function Pi=calcPiGHK(y,Z,c,GHKdraws)
        
        [T,J]=size(Z);

        L1 = tril(ones(J-1),0);
        L1(~~L1)=[1 c]; 
        L=[zeros(J-1,1) L1];
        L=[zeros(1,J); L];
        Omega=L*L';
        
        
        for t=1:T %do not parfor the GHK unless you check that the RNG is the same for each t on each call
            j=y(t);

            temp=eye(J-1); 
            Mi=[temp(:,1:j-1) -1*ones(J-1,1) temp(:,j:J-1)];

            Vi=Mi*Z(t,:)';
            %Vi=zeros(JJ-1,1);

            Li=chol(Mi*Omega(1:J,1:J)*Mi')';


            Pi(t)=ghktrain(Vi,Li,GHKdraws);   %Sim Pi 

        end
end

function Pi=calcPiLinear(y,v)
       
        logPi =  -.5*log(2*pi) - log(sigma) - .5*sigma^(-2) *(y-v(:,1))^2; %contribution of eah observation
        Pi=exp(logPi);
end