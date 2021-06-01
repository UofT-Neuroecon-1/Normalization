function out=MLestimation(dataIn,par0,opts)
    %this function sets up the ML estimation:
    %-checks data conditions
    %-sets up the dataset (pooled or hierarchical)
    %-calls optimization routine
    %
    %The LL function is also included below as a nested function so that the LL
    %has full scope (i.e. variabes do not have to be passed to it locally)
    
    %Some notes:
    %-the base altenative, i0, matters only for the beta par vector. The chol
    %decomp of the cov matrix is still parameterized in terms of the 1st item.
    
    %-If J is not constant on each trial, cannot let the covariance matrix
    %be free. It must be independent.
    

   

    
    model=opts.Models{1};
    
    options.MaxFunEvals = 2000;
   
    %% Is data unbalanced?
    T1=length(dataIn(1).y);  
    opts.balanced=1;
    for s=1:numel(dataIn)       
        if length(dataIn(s).y)~=T1
            opts.balanced=0;
            warning('Dataset is unbalanced')
        end
        
        if isrow(dataIn(s).J)
           dataIn(s).J=dataIn(s).J'; 
        end
        
    end

    %% Pooled or Clustered or Subject or Random Effect
    
    if ~isempty(opts.Hier) 
        data=dataIn;
        clear dataIn
               
    else %it is not hierarchical, so pool it
        data.X = {};
        data.Z = {};
        data.W = {};
        data.y = [];
        data.J = [];
        data.K = [];
        data.cluster=[];
        for s=1:numel(dataIn)
            data.X = [data.X,dataIn(s).X];
            if isfield(dataIn,'Z')
                data.Z = [data.Z,dataIn(s).Z];
            end
            if isfield(dataIn,'W')
                data.W = [data.W,dataIn(s).W];
            end
            data.y = [data.y;dataIn(s).y];
            data.J = [data.J;dataIn(s).J];
            data.K = [data.K;dataIn(s).K'];
            data.cluster=[data.cluster; s*ones(length(dataIn(s).X),1)];
        end
        clear dataIn
    end
    

    
    %% Choice set Properties
    opts.setsize=0; %Initially assume that choice set size doesn't change size
    for s=1:numel(data)
        T=length(data(s).y);
        Jm=max(data(s).J); %assumes each subject sees the biggest set. If not, must set opts.setsize to be max over ALL subjects.
        if ~all(data(s).J==Jm)
            opts.setsize=1; %set size is changing
        end
        
        if ~all(sort(data(s).X{1},'descend')==data(s).X{1})
            warning('Choice alternatives not ordered.')
        end
        
        if ~strcmp(opts.Prob,'Linear')
            % Pre-calculate Image Matrices for Chosen Alternative
            temp=eye(Jm-1); 
            for i=1:Jm
                M{i}=[temp(:,1:i-1) -1*ones(Jm-1,1) temp(:,i:Jm-1)];
            end
            for t=1:T
                %data(s).Mi{t}=M{data(s).y(t)}(1:data(s).J(t)-1,1:data(s).J(t));
                data(s).Mi(:,:,t)=M{data(s).y(t)}(1:data(s).J(t)-1,1:data(s).J(t));
            end  
        end
    end
     
   alg='fmincon';
   opts.objfun=@LLfun; 

%% Set Restrictions 
opts=setRestrictions(model,Jm,opts);
LB=opts.LB;
UB=opts.UB;
    
 %% Set starting values 

        
if length(par0)==sum(LB~=UB)
    theta0=par0;
else
    error('Number of initial points does not match number of parameters specified.');
end


disp('Initial Values:')    
disp(theta0)
%disp(strcat(char([opts.toEst';[opts.Hier',opts.HierDist]]),': ',strjust(num2str(num2str(theta0')),'right')))
    
if opts.getP %just getting Choice Probs, call LL and exit now
    [out.LL, out.P]=LLfun(theta0);
    out
return
end
   
 %% Start Estimation   
    disp('Start estimation');
    fprintf('# of Clusters: %d \n', max(opts.cluster));
    fprintf('# of Observations per Cluster: %d \n', T);
    
    fprintf('Size of initial sample for starting values: %d \n',opts.numInit);
    opts.Models
    opts.Prob
    opts.WithinSubject
    tic;
    
    if strcmp(alg,'optimize') 
        disp('Using Simplex Method');

        %Use optimize toolbox
        %options.Display='plot';
        options.MaxIter=10000;
        options=optimset(options,'Display',opts.Display);


        time=('00:10:00');
        getHess=0;
    elseif strcmp(alg,'fminsearch')
         disp('Using Simplex Method');
         options=optimset(options,'Display',opts.Display);

          
    elseif strcmp(alg,'fminsearchbnd')
         disp('Using Bounded Simplex Method');
         options=optimset(options,'Display',opts.Display);
         
         optfun=@fminsearchbnd;

        
    elseif strcmp(alg,'fminsearchcon')
         disp('Using (Constrained) Simplex Method');
         options=optimset(options,'Display',opts.Display,'OutputFcn',@outfun);

    elseif strcmp(alg,'fminunc')
         disp('Using Unconstrained Quasi-Newton');
         options=optimset(options,'Largescale','off','Display',opts.Display);

         getHess=1;

    elseif strcmp(alg,'fmincon')
         disp('Using Constrained Quasi-Newton');
         options=optimset(options,'Algorithm','interior-point','Display',opts.Display,'OutputFcn',@outfun);
        optfun=@fmincon;

         getHess=1;

    end

    disp(sprintf('TolX: %g   TolFun: %g',options.TolX,options.TolFun));

%     if isempty(theta0)
%         -LLfun(LB)
%     end
        
%     if size(theta0,1)==1 %start estimation at theta0 using gradiant method
%         %[thetah, maxLL, exitflags]=feval(optfun,opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
%         %[thetah, maxLL, exitflags]=fminsearchbnd(opts.objfun,theta0,LB(LB~=UB),UB(LB~=UB),options);
%         %[thetah, maxLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,thetah,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
%         
%         [thetah, nLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
%         i=1;
% 
%         maxLL=-nLL;
%         parh=thetah(i,:);
% 
%         fprintf('Value of the log-likelihood function at convergence: %9.4f \n',maxLL(i));
%         exitflag=exitflags(i);
%         disp(['Estimation took ' num2str(toc./60) ' minutes.']);
%         disp('Estimates:');
%         names=[opts.toEst opts.toRestr];
%         disp(strcat(char([opts.toEst';opts.Hier';opts.toRestr']),': ',strjust(num2str(num2str([parh'; opts.thetaR'])),'right')))
%     else
%         
        
        %theta0=randn(1, sum(LB~=UB));
        
%         [thetah, nLL, exitflags, ~]=rmsearch(opts.objfun,'fminsearchbnd',theta0,LB(LB~=UB),UB(LB~=UB),'options',options,'InitialSample',opts.numInit);
%             [~,i]=max(-nLL);
% 
%         fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-nLL(i));
%         exitflag=exitflags(i);
%         disp(['Estimation took ' num2str(toc./60) ' minutes.']);
%         disp('Estimates:');
%         disp()
%         disp(thetah);
%         save 'mpnormEst1stStage.mat' 'thetah' 'nLL' 'exitflags'

        for n=1:size(theta0,1)
           [thetah(n,:), nLL(n), exitflags(n)]=fminsearchbnd(opts.objfun,theta0(n,:),LB(LB~=UB),UB(LB~=UB),options);  
            
           [thetah(n,:), nLL(n), exitflags(n),~,~,grad(n,:),hess]=fmincon(opts.objfun,thetah(n,:),[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        end
        [maxLL,i]=max(-nLL);
        parh=thetah(i,:);
        
        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',maxLL);
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        names=[opts.toEst opts.toRestr];
        disp(strcat(char([opts.toEst';opts.Hier';opts.toRestr']),': ',strjust(num2str(num2str([parh'; opts.thetaR'])),'right')))

                       

    %Get likelihood at estimates for Vuong test	
    [~,P]=opts.objfun(thetah(i,:));
    
    %Get gradient for convergence check 
    %Or just report gradient from optimizer
    %grad=gradest(@(x) -LLfun(x),thetah(i,:));    
    
    checkConvergence;
    
    save 'mpestimates.mat' 'data' 'grad' 'Jm'

    disp('Estimates saved to disk');
    
    if opts.ses==1
        disp('Calculating finite-difference hessian and taking inverse for standard errors...');
        H=hessian(@(x) -LLfun(x),thetah(i,:));
        covh=inv(-H);
        ses=sqrt(diag(covh))';

    else
        H=[];
        covh=[];
        ses=zeros(1,length(LB));
    end
    disp('Standard Errors:');
    disp(strcat([char(opts.toEst');char(opts.Hier')],': ',strjust(num2str(num2str(ses')),'left')))
    save 'mpestimates.mat' 'data' 'H' 'grad' 'Jm' 'ses'
   
    out.exitflag=exitflag;
    out.theta0=theta0;
    out.max=i;
    out.parh=parh;
    out.LL=maxLL;
    out.P=P;
    out.se=ses;
    out.cov=covh;
    out.LB=LB;
    out.UB=UB;
    out.model=model;
    out.Prob=opts.Prob;
    out.toEst=opts.toEst;
        
    %%%%%% Nested Functions %%%%%%%   
    function [nLL, logPi]=LLfun(theta)
        
        if ~isempty(opts.Hier)
                      
            %par=opts.LB(1:6);
            
            %toEstimate=opts.LB(1:6)~=opts.UB(1:6);
            
            %par(toEstimate)=theta(1:end-sum(opts.Hier));
            %par=theta(1:opts.toEst);
            
            
            rng(1,'twister') % set random number seed
            R=opts.R; %200
            S=opts.cluster;

            for k=1:length(opts.Hier)
                par=opts.Hier(k);
                
                %p = haltonset(S,'Skip',1e3,'Leap',1e2);  %Halton Sequence
                HierParams=theta(strcmp(opts.Hier(k),[opts.toEst,opts.Hier]));
    
                %draws(:,:,k)=gaminv(net(p,R),HierParams(1),HierParams(2)); %need to draw independently for each subject for consistency and asymptotic normality
                if strcmp(opts.HierDist(k),'gamma')
                    draws(:,k,:)=gamrnd(HierParams(1),HierParams(2),S,1,R);
                elseif strcmp(opts.HierDist(k),'normal')
                    draws(:,k,:)=normrnd(HierParams(1),HierParams(2),S,1,R);
                end
            
            end

            
            
            %draws=reshape(logninv(net(p,R*S),par(opts.Hier==1),theta(end-sum(opts.Hier)+1:end)),R,S);
            %if opts.balanced
            
                tic
                %Pi=zeros(length(data(s).X),R,S);
                logPi=cell(S,1);
                
                parfor subj=1:S 
    %                 for r=1:R
    %                     particle=struct();
    %                     par=theta(1:length(opts.toEst));                   
    %                     par(strcmp(opts.Hier,opts.toEst))=squeeze(draws(r,subj,:));
    %                     particle.theta=par;
    %                     Pi(:,subj,r)=ProbaChoice(data(subj), particle, opts );            
    %                 end

                        particle=struct();
                        
                        par=repmat(theta(1:length(opts.toEst)),1,1,R);     %just initialize              
                        par(:,strcmp(opts.Hier,opts.toEst),:)=draws(subj,:,:);
                        particle.theta=par;
                       
                        logPi{subj}=ProbaChoice(data(subj), particle, opts );
                        %Pi(:,:,subj)=ProbaChoice(data(subj), particle, opts );            

                end
                toc
                %nLL=-sum(log(mean(prod(Pi),2))); Too many zeros, line below corrects numerical issues
                
                nLLs=@(x) (log(mean(exp(T*log(3)+sum(x)),2)) - T*log(3));            
                nLL=-sum( cellfun(nLLs,logPi) );
                
                %nLL=-sum( log(mean(exp(T*log(3)+sum(log(Pi))),2)) - T*log(3) );
%             else
%                 error('Dataset isnt balanced, cannot compute LL (yet)');
%             end
%             
        else
            particle.theta=theta;
            logPi=ProbaChoice(data, particle, opts );
            nLL=-sum(logPi);
        end
    end
            
    function checkConvergence
        disp(' ');

        if exitflag == 1
            disp('Convergence achieved.');

        elseif exitflag == 2
            disp('Convergence achieved by criterion based on change in parameters.');
            if size(options.TolX,1)>0
                disp(['Parameters changed less than PARAMTOL= ' num2str(options.TolX)]);
            else
                disp('Parameters changed less than PARAMTOL=0.000001, set by default.');
            end
            disp('You might want to check whether this is actually convergence.');
            disp('The gradient vector is');
            disp(grad)
        elseif exitflag == 3
            disp('Convergence achieved by criterion based on change in log-likelihood value.');
            if size(options.TolFun,1)>0
                disp(['Log-likelihood value changed less than LLTOL= ' num2str(options.TolFun)]);
            else
                disp('Log-likelihood changed less than LLTOL=0.000001, set by default.');
            end
            disp('You might want to check whether this is actually convergence.');
            disp('The gradient vector is');
            disp(grad)

        elseif exitflag == 5
            disp('Predicted decrease in the objective function was less than the TolFun tolerance.');
            disp('Convergence likely achieved, but the first-order optimality measure was above the function tolerance.');

        else
            disp('Convergence not achieved.');
            disp('Results are not printed because no convergence.');
            return
        end
    end

    function stop = outfun(x,~,state) %display parameter values during estimation
        stop = false; 
        if state == 'iter' 
            disp([' x = ',num2str(x)]);
            save('iteroutput')
        end 
    end
end