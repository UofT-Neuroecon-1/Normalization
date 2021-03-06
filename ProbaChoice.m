function [ Pi ] = ProbaChoice( data, particle,opts )
    % data.X: J x K matrix of choice set
    % data.y: Tx1 of choices
    % model: the model to use
    % (returns) Pi : Jx1 vector of choice probabilities
    %F: function handle for model
    

    J = data.J;
    K = data.K;
    
    F=eval(['@' opts.modelF]);
        
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
    
    Pi=F(data.X,data.Z,data.W,particle.theta); %Get Probs
 
    Pi(Pi==0)=realmin;  % make sure that the probability of a choice is not 0


    function Pi=DN(X,~,~,theta)
        
        par=[theta repmat(opts.thetaR,size(theta,1),1)];
        names=[opts.toEst opts.toRestr];
        
        s = par(:,strcmp(names, 'sigma'));
        w = par(:,startsWith(names, 'omega')); % all weights
        a = par(:,strcmp(names, 'alpha'));
        wx= par(:,strcmp(names, 'wx')); %additional own weight
        
        f = @(x) (x.^(a'));
        %denom=@(x) (s + w*sum(x) );
        
        %
        if any(strcmp(opts.toEst, 'beta'))
            b = par(strcmp(names, 'beta'));
            if b<0
                b=0.0000000001; %hack so that s.e. calculation (hessian()) can send negative b. In estimation, beta is restricted to be >0 so doesn't matter.
            end
            weights=w*ones(3,3);
            ind=eye(3);
            if any(strcmp(opts.toEst, 'wx'))
                weights(logical(ind(:)))=wx;
            end
            
            denom=@(x) (s + vecnorm(weights.*x',b,2) ); %this line includes the weight in the norm, so (w*x)^b.
            %denom=@(x) (s + sum(weights.*(x.^b'),2).^(1/b) ); %messed up
            
            %Note that beloe wx is the "extra" weight relative w. So in
            %reporting results, add them together for w_i. Can be messed up when w<0 and taking norm (i.e. abs())
            %denom=@(x) (s + vecnorm((w + wx*eye(length(x))).*x',b,2) ); %this line includes the weight in the norm, so (w*x)^b.
            %denom=@(x) (s + sum(((w + wx*eye(length(x))).*x).^b',2).^(1/b) ); %this line does not include the abs(weight) in the norm, so w*x^b.
            %denom=@(x) (s + sum((w + wx*eye(length(x))).*(x.^b'),2).^(1/b) ); %this line does not include the abs(weight) in the norm, so w*x^b.
        elseif any(strcmp(opts.toEst, 'wx'))  
            weights=w*ones(3,3);
            ind=eye(3);
            if any(strcmp(opts.toEst, 'wx'))
                weights(logical(ind(:)))=wx;
            end
            denom=@(x) (s + vecnorm(weights.*x',1,2) ); %this line is not wrong, it just doesn't allow beta free, but allows w<0.
        else
            if size(w,2)==1
                weight=@(x) repmat(w,1,length(x));
            else
                weight=@(x) w;
            end
            denom=@(x) (s + weight(x)*x)'; %if w is a scalar, then expand to vector and take vector product
        end
        
        
        %denom=@(x) (sigma + omega*norm(x,b) );
        %vecnorm(cell2mat(X)',2)
        sumv=cellfun(denom,X,'uniformoutput',false);
        numerator=cellfun(f,X,'uniformoutput',false);
        v=cellfun(@rdivide,numerator,sumv,'uniformoutput',false);

        
        if ~strcmp(opts.Prob,'GHK')
            P=eval(['@calcPi' opts.Prob]); %Construct function handle for probability function based on label in opts.Prob
            Pi = P(data.Mi,v,data.J,opts.scale); %y is Mi
        else %use GHK

            if opts.setsize==1
                Pi=calcPiGHKC(data.Mi,v,data.J,par(find(startsWith(names, 'c'))));
            else
                Pi=calcPiGHK(data.y,cell2mat(v)',par(find(startsWith(names, 'c'))),opts.GHKdraws);
            end
        end
        
    end

    function Pi=Range(X,~,~,theta)
 
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
        
        denom=@(x) (s + w*(max(x)-w2*min(x)));
        
        %denom=@(x) (sigma + omega*norm(x,b) );
        %vecnorm(cell2mat(X)',2)
        sumv=cellfun(denom,X,'uniformoutput',false);
        numerator=cellfun(f,X,'uniformoutput',false);
        v=cellfun(@rdivide,numerator,sumv,'uniformoutput',false);
        
        P=eval(['@calcPi' opts.Prob]);
        Pi = P(data.Mi,v,data.J); %y is Mi
            
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
        wx= par(find(strcmp(names, 'webb'))); %additional own weight
        
       
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
        
        P=eval(['@calcPi' opts.Prob]);
        Pi = P(data.Mi,eu,data.J); %y is Mi
    end
end

%Get Choice Probs given the model above

function Pi=calcPiProbit(Mi,v,J,scale)

    T=size(v,2);
    
    
    if length(v{1})==2 %Binary Choice
        if all(J==J(1)) % all choice sets same size, then all trials together (much faster)   
            vi=cellfun(@mtimes, Mi, v)';
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
end

function Pi=calcPiLogit(Mi,v,~,scale)
        vi=cellfun(@mtimes, Mi, v, 'UniformOutput', false);
        %cellfun(@exp,vi','UniformOutput', false)
        sum_exp_v = cellfun(@(x) sum(exp(x./scale),1),vi,'UniformOutput', false);
        Pi = 1./(1+cell2mat(sum_exp_v'));
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