    function v=DN(data,par,names,opts)
       
        
        s = par(:,strcmp(names, 'sigma'),:);
        w = par(:,startsWith(names, 'omega'),:); % all weights
        a = par(:,strcmp(names, 'alpha'),:);
        wx= par(:,strcmp(names, 'wx'),:); %additional own weight
        v0= par(:,strcmp(names, 'v0'),:);
        
        %f = @(x) (x(sum(opts.toNorm,2)>1).^(a'));
        
        %
        if any(strcmp(opts.toEst, 'beta'))
            b = par(strcmp(names, 'beta'));
            if b<0
                b=0.0000000001; %hack so that s.e. calculation (hessian()) can send negative b. In estimation, beta is restricted to be >0 so doesn't matter.
            end
            weights=w*opts.toNorm;
            ind=eye(3); error('Need to adjust index to match number of alternatives and conform with opts.toNorm')
            if any(strcmp(opts.toEst, 'wx'))
                weights(logical(ind(:)))=wx;
            end
            
            %denom=@(x) (s + vecnorm(weights.*x',b,2) ); %this line includes the weight in the norm, so (w*x)^b.
            
            V=@(x) ( (x(sum(opts.toNorm,2)>1,:).^(a')) ./ (s + vecnorm(weights.*x',b,2) ));
            %denom=@(x) (s + sum(weights.*(x.^b'),2).^(1/b) ); %messed up
            
            %Note that beloe wx is the "extra" weight relative w. So in
            %reporting results, add them together for w_i. Can be messed up when w<0 and taking norm (i.e. abs())
            %denom=@(x) (s + vecnorm((w + wx*eye(length(x))).*x',b,2) ); %this line includes the weight in the norm, so (w*x)^b.
            %denom=@(x) (s + sum(((w + wx*eye(length(x))).*x).^b',2).^(1/b) ); %this line does not include the abs(weight) in the norm, so w*x^b.
            %denom=@(x) (s + sum((w + wx*eye(length(x))).*(x.^b'),2).^(1/b) ); %this line does not include the abs(weight) in the norm, so w*x^b.
        elseif any(strcmp(opts.toEst, 'wx'))  
            weights=w*opts.toNorm;
            ind=eye(3); error('Need to adjust index to match number of alternatives and conform with opts.toNorm')
            if any(strcmp(opts.toEst, 'wx'))
                weights(logical(ind(:)))=wx;
            end
            V=@(x) ( x(sum(opts.toNorm,2)>1,:).^(a')) ./ (s + vecnorm(weights(sum(opts.toNorm,2)>1,:).*x',1,2) );
            %denom=@(x) (s + vecnorm(weights(sum(opts.toNorm,2)>1,:).*x',1,2) ); %this line is not wrong, it just doesn't allow beta free, but allows w<0.
        else
            if size(w,2)==1
                weights=reshape(w,1,1,length(w)).*opts.toNorm;
            else
                weights=w; warning('Check dimension of vector multiplication');
            end
            V=@(x) ( (x(sum(opts.toNorm,2)>1,:).^(a)) ./ (reshape(s,1,1,length(s)) + squeeze(pagemtimes(weights(sum(opts.toNorm,2)>1,:,:),x))) );
            %denom=@(x) (s' + squeeze(pagemtimes(weights(sum(opts.toNorm,2)>1,:,:),x))); %if w is a scalar, then expand to vector and take vector product
        end
        
       
%         sumv=cellfun(denom,X,'uniformoutput',false);
%         numerator=cellfun(f,X,'uniformoutput',false);
%         v=cellfun(@rdivide,numerator,sumv,'uniformoutput',false);

        %v=cellfun(V,X,'uniformoutput',false);
        v=v0+V(data.X);
    end