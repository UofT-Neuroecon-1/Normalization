function out=MLestimation(dataIn,par0,opts)

    %the base altenative, i0, matters only for the beta par vector. The chol
    %decomp of the cov matrix is still parameterized in terms of the 1st item.
    
    %If J is not constant on each trial, cannot let the covariance matrix
    %be free. It must be independent.
    

   

    
    model=opts.Models{1};
    
    options.MaxFunEvals = 2000;
   
    

    %% Pooled or Clustered or Subject or Random Effect
    
    if any(opts.Hier) 
        data=dataIn;
        clear dataIn
    else %it is not hierarchical
        data.X = {};
        data.y = [];
        data.J = [];
        data.K = [];
        data.cluster=[];
        for s=1:numel(dataIn)
            data.X = [data.X,dataIn(s).X];
            data.y = [data.y;dataIn(s).y];
            data.J = [data.J;dataIn(s).J'];
            data.K = [data.K;dataIn(s).K'];
            data.cluster=[data.cluster; s*ones(length(dataIn(s).X),1)];
        end
        clear dataIn
    end
    
        %% Choice set Properties
    for s=1:numel(data)
        T=size(data(s).X,2);

        Jm=max(data(s).J); %assumes each subject sees the biggest set, so is set here. If not, must set opts.setsize to be max over ALL subjects.
        if ~all(data(s).J==Jm)
            opts.setsize=1; %set size is changing
        end
        
        % Pre-calculate Image Matrices for Chosen Alternative
        temp=eye(Jm-1); 
        for i=1:Jm
            M{i}=[temp(:,1:i-1) -1*ones(Jm-1,1) temp(:,i:Jm-1)];
        end
        for t=1:T
            data(s).Mi{t}=M{data(s).y(t)}(1:data(s).J(t)-1,1:data(s).J(t)); 
        end      
    end

    
if isfield(data,'W')
    W=data.W;
else W=[];
end
    

    
%     if indep==0
%         alg='fminunc';  
%     else %and homoskedastic
%         %alg='fmincon';
%         alg='fminsearchbnd';
%     end
   % alg='fminsearchbnd';
   alg='fmincon';
   opts.objfun=@LLfun; 

%% Set Restrictions 
opts=setRestrictions(model,Jm,opts);
LB=opts.LB;
UB=opts.UB;
    
 %% Set starting values 
    if isempty(par0)
        disp('No Initial parameters specified, starting point is random');

        theta0=randn(1, sum(LB~=UB));
    else
        theta0=par0;
    end
disp('Initial Values:')    
disp(theta0)
    
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
    tic;
    

    

    if strcmp(alg,'optimize') 
        display('Using Simplex Method');

        %Use optimize toolbox
        %options.Display='plot';
        options.MaxIter=10000;
        options=optimset(options,'Display',opts.Display);

%         optfun=@optimize;
%         varin={[],[],[],[],[],[],[],[],options,alg};

         time=('00:10:00');
        getHess=0;
    elseif strcmp(alg,'fminsearch')
         display('Using Simplex Method');
         options=optimset(options,'Display',opts.Display);
 %        optfun=@fminsearch;
 %        varin={options};
          
    elseif strcmp(alg,'fminsearchbnd')
         display('Using Bounded Simplex Method');
         options=optimset(options,'Display',opts.Display);
         
          optfun=@fminsearchbnd;
 %       varin={LB(LB~=UB),UB(LB~=UB),options};
        
    elseif strcmp(alg,'fminsearchcon')
         display('Using (Constrained) Simplex Method');
         options=optimset(options,'Display',opts.Display,'OutputFcn',@outfun);
%          optfun=@fminsearchcon;
%        varin={LB(LB~=UB),UB(LB~=UB),LIN(LB~=UB),0,[],options};

    elseif strcmp(alg,'fminunc')
         display('Using Unconstrained Quasi-Newton');
         options=optimset(options,'Largescale','off','Display',opts.Display);
 %        optfun=@fminunc;
 %        varin={options};

         getHess=1;

    elseif strcmp(alg,'fmincon')
         display('Using Constrained Quasi-Newton');
         options=optimset(options,'Algorithm','interior-point','Display',opts.Display,'OutputFcn',@outfun);
        optfun=@fmincon;
%          varin={[],[],LIN(LB~=UB),0,LB(LB~=UB),UB(LB~=UB),[],options};
%        varin={[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options};

         getHess=1;

    end

    display(sprintf('TolX: %g   TolFun: %g',options.TolX,options.TolFun));

%     if isempty(theta0)
%         -LLfun(LB)
%     end

    if opts.numInit==1
        %[thetah, maxLL, exitflags]=feval(optfun,opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        %[thetah, maxLL, exitflags]=fminsearchbnd(opts.objfun,theta0,LB(LB~=UB),UB(LB~=UB),options);
        %[thetah, maxLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,thetah,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        
        [thetah, maxLL, exitflags,~,~,grad,hess]=fmincon(opts.objfun,theta0,[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        i=1;

        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-maxLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
    else
        [thetah, nLL, exitflags, xstart]=rmsearch(opts.objfun,'fminsearchbnd',theta0,LB(LB~=UB),UB(LB~=UB),'options',options,'InitialSample',opts.numInit);
            [maxLL,i]=max(-nLL);

        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-nLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
        save 'mpnormEst1stStage.mat' 'thetah' 'nLL' 'exitflags'

        for n=1:opts.numInit
            [thetah(n,:), nLL(n), exitflags(n),~,~,grad(n,:),hess]=fmincon(opts.objfun,thetah(n,:),[],[],[],[],LB(LB~=UB),UB(LB~=UB),[],options);
        end
        [maxLL,i]=max(-nLL);
 
        fprintf('Value of the log-likelihood function at convergence: %9.4f \n',-nLL(i));
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        disp(thetah);
        save 'mpnormEst2ndStage.mat' 'thetah' 'nLL' 'exitflags' 'grad'
    end
                       
%         if strcmp(alg,'fminsearch') || strcmp(alg,'fminsearchcon') || strcmp(alg,'fminsearchbnd') || strcmp(alg,'optimize')
%             [thetah(s,:), nLL(s,1), exitflags(s,1)]=optfun(objfun,theta0(s,:),varin{:});
%             grad=[]
%             hess=[];
%         else
%                 if strcmp(alg,'interior-point')
%                     
%                     [xfinal,ffinal,exitflag,xstart] = rmsearch(fun,optname,x0,LB(LB~=UB),UB(LB~=UB),options)
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),xstart]=rmsearch(objfun,optfun,theta0(s,:),LB(LB~=UB),UB(LB~=UB),options);
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),~,~,grad,H(:,:,s)]=optfun(objfun,theta0(s,:),varin{:});
%                 else
%                     [thetah(s,:), nLL(s,1), exitflags(s,1),~,grad,H(:,:,s)]=optfun(objfun,theta0(s,:),varin{:});
%                 end
% 
%         end


    %matlabpool close
   

    %Get likelihood at estimates for Vuong test	
    [~,P]=opts.objfun(thetah(i,:));
    
    %Get gradient for convergence check 
    %Or just report gradient from optimizer
    %grad=gradest(@(x) -LLfun(x),thetah(i,:));    
    
    checkConvergence;


    % if PREDICT == 1;
    %     disp(' ');
    %     disp('Predict shares at estimated coefficients.');
    %     probs=pred(paramhat);
    % end;
    i0=opts.i0;      
    save 'mpestimates.mat' 'data' 'grad' 'Jm' 'i0'

    display('Estimates saved to disk');
    
    if opts.ses==1
    disp('Calculating finite-difference hessian and taking inverse for standard errors.');
    H=hessian(@(x) -LLfun(x),thetah(i,:));
    covh=inv(-H);
    se=sqrt(diag(covh))';

    ses(LB~=UB)=se;
    ses(LB==UB)=0;
    else
        H=0;
        covh=inv(-H);
       ses=zeros(1,length(LB));
    end
    disp('Standard Errors:');
    disp(ses);
    save 'mpestimates.mat' 'data' 'H' 'grad' 'Jm' 'ses' 'i0'
    
   parh=LB; 
   parh(LB~=UB)=thetah(i,:);
   
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
        
    %%%%%% Nested Functions %%%%%%%   
    function [nLL, Pi]=LLfun(theta)
        
        if any(opts.Hier)
           
            R=200;
            S=opts.cluster;
            
            par=opts.LB(1:6);
            
            toEstimate=opts.LB(1:6)~=opts.UB(1:6);
            
            par(toEstimate)=theta(1:end-sum(opts.Hier));

            rng(1,'twister') % set random number seed
            
            %draws=gamrnd(par(opts.Hier==1),theta(end-sum(opts.Hier)+1:end),R,1);
            
            p = haltonset(1,'Skip',1e3,'Leap',1e2);  %Halton Sequence
            draws=reshape(gaminv(net(p,R*S),par(opts.Hier==1),theta(end-sum(opts.Hier)+1:end)),R,S); %get draws
            %draws=reshape(logninv(net(p,R*S),par(opts.Hier==1),theta(end-sum(opts.Hier)+1:end)),R,S);
            
            Pi=zeros(length(data(s).X),S,R);

            parfor s=1:S
            for r=1:R
                particle=struct();
                temp2=par;
                temp2(opts.Hier)=draws(r,s);
                particle.theta=temp2(toEstimate);
                Pi(:,s,r)=ProbaChoice(data(s), particle, opts );            
            end
            end

            %nLL=-sum(log(mean(prod(Pi),2))); Too many zeros, line below corrects numerical issues
            nLL=-sum(log(mean(exp(T*log(3)+sum(log(Pi))),3)) - T*log(3));
            
        else
            particle.theta=theta;
            Pi=ProbaChoice(data, particle, opts );
            nLL=-sum(log(Pi));
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

    function stop = outfun(x,~,state) 
        stop = false; 
        if state == 'iter' 
            disp([' x = ',num2str(x)]); 
        end 
    end
end