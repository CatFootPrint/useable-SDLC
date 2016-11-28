function [D,Y,Out] = dl_proxi_advanced(X,K,lambda1,lambda2,opts)
%learn a dictionary D from
%---- min_{D,Y} 0.5*||X-D*Y||_F^2 + mu*||Y||_1 --------------%
%---- subject to norm(D(:,j)) <= 1 for all j ----------------%

[m,p] = size(X); 

if isfield(opts,'tol')       tol = opts.tol;      else tol = 1e-4;     end
if isfield(opts,'maxit')     maxit = opts.maxit;  else maxit = 500;    end
if isfield(opts,'maxT')      maxT = opts.maxT;    else maxT = 1e3;     end
% yType = 0: no nonnegativity constraint
% yType = 1: Y is nonnegative
if isfield(opts,'yType')     yType = opts.yType;  else yType = 0;      end
if isfield(opts,'D0')        D0 = opts.D0;        else D0 = randn(m,K);end
if isfield(opts,'Y0')        Y0 = opts.Y0;        else Y0 = randn(K,p);end

% normalize D0
for j=1:K
    D0(:,j) = D0(:,j)/norm(D0(:,j));
end

nrmX = norm(X,'fro');
D = D0; Dm = D0; Dt = D'; Dsq = Dt*D;
Y = Y0; Ym = Y0; Yt = Y'; Ysq = Y*Yt;
t0 = 1; t = 1; Ly = 1; Ld = 1; 
rw = 0.9999;

obj0 = 0.5*norm(D0*Y0-X,'fro')^2+lambda1*sum(sum(abs(Y0)));
obj = obj0;

Out.redoN = 0; nstall = 0; start_time = tic;
fprintf('Iteration of advanced SDLC:     ');

for k = 1:maxit
    fprintf('\b\b\b\b\b%5i',k);
    Out.iteration_times=k;
    % update Y
    Ly0 = Ly;       Ly = norm(Dsq); 
    GradY = Dt*(D*Ym-X)+(lambda2)*Ym;
    Phi = Ym-GradY/Ly;
    Y=Phi-1/4/Ly*sign(Phi-lambda1/2/size(D0,2)/Ly*ones(size(Phi)))-1/4/Ly*sign(Phi+lambda1/2/size(D0,2)/Ly*ones(size(Phi)));
    % update D
    Ld0 = Ld;       
%     Ld = norm(Ysq);
    for i=1:size(D0,2)
        ai=transpose(Y(i,:));
%         di=Dm(:,i);
        Ld = norm(ai); 
%         Ei=X-Dm*Y+di*ai';
        GradD = (X-Dm*Y)*ai;
        D(:,i) = D(:,i)-GradD/Ld;
    end
        nrmD = sqrt(sum(D.^2,1));
        id = nrmD>1; s = sum(id);
        if s>0     
            D(:,id) = D(:,id)*spdiags(nrmD(id)'.^(-1),0,s,s);
        end
    Dt = D'; Dsq0 = Dsq; Dsq = D*Dt;
    
    res = norm(D*Y-X,'fro');
    obj0 = obj; 
    switch yType
        case 0
            obj = res^2+lambda1*sum(abs(Y(:)))+(lambda2)*sum((abs(Y(:))).^2);
        case 1
            obj = res^2+lambda1*sum(sum(Y))+(lambda2)*sum((abs(Y)).^2);
    end
    
    
    if obj>obj0
        Out.redoN = Out.redoN+1;
        
       % update Y
        Ly = norm(Dsq0);
        GradY = D0'*(D0*Y0-X);
        GradY = D0'*(D0*Y0-X)+(lambda2)*Y0;
        Y = Y0-GradY/Ly;




%         GradY = D0'*(D*Y0-X)+(1-mu)*Y0;
%     Phi = Y0-GradY/Ly;
%     Y=Phi-1/4/Ly*sign(Phi-mu/2/size(D0,2)/Ly*ones(size(Phi)))-1/4/Ly*sign(Phi+mu/2/size(D0,2)/Ly*ones(size(Phi)));
        
        switch yType
            case 0
                Y = sign(Y).*max(abs(Y)-lambda1/Ly,0);
            case 1
                Y = max(Y-lambda1/Ly,0);
        end  
        
        Yt = Y';
        Ysq = Y*Yt;
        
            
        
        % update D
        Ld = norm(Ysq);
        GradD = D0*Ysq-X*Yt;
        D = D0-GradD/Ld;
        nrmD = sqrt(sum(D.^2,1));
        id = nrmD>1; s = sum(id);
        if s>0     
            D(:,id) = D(:,id)*spdiags(nrmD(id)'.^(-1),0,s,s);
        end
        Dt = D'; Dsq = Dt*D;
        
        res = norm(D*Y-X,'fro');
        obj0 = obj; 
        switch yType
            case 0
%                 obj = 0.5*res^2+mu*sum(abs(Y(:)));
obj = res^2+lambda1*sum(abs(Y(:)))+(lambda2)*sum((abs(Y(:))).^2);
            case 1
%                 obj = 0.5*res^2+mu*sum(sum(Y));
obj = res^2+lambda1*sum(sum(Y))+(lambda2)*sum((abs(Y)).^2);
        end
        Out.D{k} = D;
    end
    % do extrapolation
    t = (1+sqrt(1+4*t0^2))/2;
    w = (t0-1)/t; % extrapolation weight
    wY = min([w,rw*sqrt(Ly0/Ly)]);
    wD = min([w,rw*sqrt(Ld0/Ld)]);
    Ym = Y+wY*(Y-Y0);   Dm = D+wD*(D-D0);
    t0 = t; Y0 = Y; D0 = D;
    
    % --- diagnostics, reporting, stopping checks ---
    relerr1 = abs(obj-obj0)/(obj0+1);    relerr2 = res/nrmX;
    
    % reporting
    Out.hist_obj(k) = obj;
    Out.hist_rel(1,k) = relerr1;
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion
    crit = relerr1<tol;
    if crit; nstall = nstall+1; else nstall = 0; end
    if nstall>=3 || relerr2<tol 
        break; 
    end
    if toc(start_time)>maxT; break; end;    
end
fprintf('\n'); 
Out.iter = k;
