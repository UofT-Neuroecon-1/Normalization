function out=LSestimation(dataIn,par0,opts)
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
    
    if iscell(data.X)
        data.X=cell2mat(data.X);
        data.Z=cell2mat(data.Z);
        data.W=cell2mat(data.W);
    end


    
    %% Choice set Properties
    opts.setsize=0; %Initially assume that choice set size doesn't change size
    for s=1:numel(data)
        T=length(data(s).y);
        Jm=max(data(s).J); %assumes each subject sees the biggest set. If not, must set opts.setsize to be max over ALL subjects.
        if ~all(data(s).J==Jm)
            opts.setsize=1; %set size is changing
        end
        
        if ~all(sort(data(s).X(1,:),'descend')==data(s).X(1,:))
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
     
   alg='lsqnonlin';
   opts.objfun=@ErrorFun; 

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
    
    if strcmp(alg,'lsqnonlin') 
        disp('Using Non-Linear Least Squares');

        %Use optimize toolbox
        %options.Display='plot';
        options.MaxIter=10000;
        options=optimset(options,'Display',opts.Display);


        time=('00:10:00');
        getHess=0;
    end

    disp(sprintf('TolX: %g   TolFun: %g',options.TolX,options.TolFun));
    
        for n=1:size(theta0,1)
           [thetah(n,:), SSE(n),Rs(:,n), exitflags(n)]=lsqnonlin(opts.objfun,theta0(n,:),LB(LB~=UB),UB(LB~=UB),options);  
            
        end
        [minSSE,i]=min(SSE);
        residuals=Rs(:,i);
        parh=thetah(i,:);
        
        fprintf('Value of the error function at convergence: %9.4f \n',minSSE);
        exitflag=exitflags(i);
        disp(['Estimation took ' num2str(toc./60) ' minutes.']);
        disp('Estimates:');
        names=[opts.toEst opts.toRestr];
        disp(strcat(char([opts.toEst';opts.Hier';opts.toRestr']),': ',strjust(num2str(num2str([parh'; opts.thetaR'])),'right')))

                       

    %Get likelihood at estimates for Vuong test	
    if ~strcmp(opts.Prob,'Linear')
        [~,P]=opts.objfun(thetah(i,:));
        
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
        
    else 
        
        ci = nlparci(parh,resisuals,'covar',sigma)

        ses=zeros(1,length(LB));

        disp('Standard Errors:');
        disp(strcat([char(opts.toEst');char(opts.Hier')],': ',strjust(num2str(num2str(ses')),'left')))
        save 'mpestimates.mat' 'data' 'Jm' 'ses'
        
        out.exitflag=exitflag;
        out.theta0=theta0;
        out.max=i;
        out.parh=parh;
        out.SSE=minSSE;
        out.se=ses;

        out.LB=LB;
        out.UB=UB;
        out.model=model;
        out.Prob=opts.Prob;
        out.toEst=opts.toEst;
    end
    
    %Get gradient for convergence check 
    %Or just report gradient from optimizer
    %grad=gradest(@(x) -LLfun(x),thetah(i,:));    
    
    
    
        
    %%%%%% Nested Functions %%%%%%%   
    function [errors]=ErrorFun(theta)
    
        par=[theta repmat(opts.thetaR,1,1,size(theta,3))];
        names=[opts.toEst opts.toRestr];

        v=DN(data, par, names, opts);
        errors=data.y - v(1,:)';
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