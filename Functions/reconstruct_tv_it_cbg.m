function [reconstruction] = reconstruct_tv_it_cbg(A,cscale, measurement, reg, NumBlocks, tik_reg)
%% Forms the scene reconstruction given a measurement image and using the A matrices.
  
   % Difference matrix
   Diff = diag(ones(1,size(A,1)),0) - diag(ones(1,size(A,1)-1),1);

   % Difference measurements
    mr = (Diff*reshape(measurement(:,:,1),[],1));
    mg = (Diff*reshape(measurement(:,:,2),[],1));
    mb = (Diff*reshape(measurement(:,:,3),[],1));
    
    dA = (Diff*A);
    AT = dA';
    ATA = AT*dA;

    % Initial estimates
    disp(' Making initial estimate...')
    
    if tik_reg>0  %Tikhinov regularization with 1/f DCT prior can be a good initialization
        weight1 = (linspace(0.1,1,NumBlocks(2)).^0.35);
        weight2 = (linspace(0.1,1,NumBlocks(1)).^0.35);
        weight1 = weight1./max(weight1);
        weight2 = weight2./max(weight2);
        D = kron(diag(weight1) *dctmtx(NumBlocks(2)),diag(weight2)*dctmtx(NumBlocks(1)));
        DtD = D'*D;
    else
        DtD = 0;
    end
    
    testr = rs(( cscale(1)^2 * ATA + tik_reg(1)^4*(DtD))\(cscale(1)*AT) * mr(:));
    testg = rs(( cscale(2)^2 * ATA + tik_reg(2)^4*(DtD))\(cscale(2)*AT) * mg(:));
    testb = rs(( cscale(3)^2 * ATA + tik_reg(3)^4*(DtD))\(cscale(3)*AT) * mb(:));
    

    %%
    % TV Parameters
    tvparams.constraints = 'none';

    % Cost function (g(x) + h(x)), gradient (dg(x)) and proximal operator
    grad = @(x,constants) rs((-constants.ATm + constants.ATA*x(:))); 
    prox = @(x,params) co_prox_tv(rs(x), params.lambda, tvparams);
    calc_H = @(x,constants) constants.lamb* sum(sum(sqrt([diff(x,1,1).^2;zeros(1,size(x,2))] +[diff(x,1,2).^2,zeros(size(x,1),1)])));
    calc_G = @(x,constants)  0.5*norm(constants.meas - constants.A*x(:)).^2;
    
    % FISTA solver params

    % Start reconstruction
    disp(' Reconstructing red channel...')
    constants.A = cscale(1)*dA;
    constants.ATm = cscale(1)*AT*mr(:);
    constants.ATA = cscale(1)^2 * ATA;
    constants.meas = mr(:);
    constants.lamb = reg(1);
    [r_opt] = fista_backtrack(grad, prox, testr, 1e-5, 1,  struct('lambda',reg(1)), calc_H, calc_G, constants);
    
    disp(' Reconstructing green channel...')
    constants.A = cscale(2)*dA;
    constants.ATm = cscale(2)*AT*mg(:);
    constants.ATA = cscale(2)^2 * ATA;
    constants.meas = mg(:);
    constants.lamb = reg(2);
    [g_opt] = fista_backtrack(grad, prox, testg, 1e-5, 1,  struct('lambda',reg(2)), calc_H, calc_G, constants);
    
    disp(' Reconstructing blue channel...')
    constants.A = cscale(3)*dA;
    constants.ATm = cscale(3)*AT*mb(:);
    constants.ATA = cscale(3)^2 * ATA;
    constants.meas = mb(:);
    constants.lamb = reg(3);
    [b_opt] = fista_backtrack(grad, prox, testb, 1e-5, 1,  struct('lambda',reg(3)), calc_H, calc_G, constants);
    

    reconstruction(:,:,1) = r_opt;%opt{1}; 
    reconstruction(:,:,2) = g_opt;%opt{2};
    reconstruction(:,:,3) = b_opt;%opt{3};
    
    
    function [image] = rs(input)
        image = reshape(input, NumBlocks(1),NumBlocks(2));
    end
   

end




function [solution] = co_prox_tv(input, lambda, params)
%% Solves prox_tv(Y) = argmin_X  1/2 ||X-Y||^2 + lambda*TV(X)
% Using method from Beck & Teboulle (2009)
% TV = Isotropic TV with reflexive boundary conditions, i.e d/dx(x>=size(input)) = 0;

% input = input matrix/image Y
% Parameters
%    params.tol         => Stopping tolerance, stops when:
%                          (error(i)-error(i-1))/error(i)<param.tol
%                          Defaults to 1e-6
%    params.constraint  => 'none' - default
%                          'non-neg'
%                          'bounds'
%    params.usegpu      => Use GPU arrays. 0 = No, 1 = Yes. Default = 0;

%% Optional inputs
if nargin==2
    params.tvtol = 1e-6;
    params.constraint = 'none';
    params.usegpu = 0;
else
    if ~isfield(params, 'tol')
        params.tvtol = 1e-6;
    end

    if ~isfield(params, 'constraint')
        params.constraint = 'none';
    end
    
    if ~isfield(params, 'usegpu')
        params.usegpu=0;
    end
end

%% Initialisation and operator function definitions
if strcmp(params.constraint,'non-neg')
    Pc = @(x) max(0,x);  %Project onto constrained set i.e non-negatvity contraints 
elseif strcmp(params.constraint,'bounds')
    Pc = @(x) max(0,min(1,x));  %Project onto constrained set [0,1]
else
    Pc = @(x) x;
end

 function [p,q] = Lt(x)
    %LT(x) = (p, q)
    %Where:   p[i,j] = x[i,j] - x[i+1,j]     q[i,j] = x[i,j] - x[i,j+1]
    %Used for [Pc(b-gamma*L(r,s))]
    
    if (params.usegpu==0)
        p = zeros(size(x));
        q = zeros(size(x));
    else
        p = zeros(size(x),'gpuArray');
        q = zeros(size(x),'gpuArray');
    end
    
    p(1:end-1,:) = -diff(x,1,1);
    q(:,1:end-1) = -diff(x,1,2);
    p(p~=p) = 0; %Nans to 0
    q(q~=q) = 0;
 end

 function [I] = L(dx,dy)
     %L(p, q)[i,j] = p[i,j] + q[i,j] - p[i-1,j] - q[i,j-1], i = 1,...,m, j = 1,...,n,

    if (params.usegpu==0)
        Ix = zeros(size(dx));
        Iy = zeros(size(dy));
    else
        Ix = zeros(size(dx),'gpuArray');
        Iy = zeros(size(dy),'gpuArray');
    end
    
    Ix(2:end,:) = diff(dx,1,1);
    Iy(:,2:end) = diff(dy,1,2);

    Ix(1,:) = dx(1,:,:);
    Iy(:,1) = dy(:,1,:);

    I = Ix + Iy;
        
 end

%% Initialise values

if (params.usegpu==0)
    b = input;
    t_old = 1;
    r = zeros(size(b));
    s = zeros(size(b));
    p_old = r;
    q_old = s;
else
    b = gpuArray(input);
    r = zeros(size(b),'gpuArray');
    s = zeros(size(b),'gpuArray');
    p_old = zeros(size(b),'gpuArray');
    q_old = zeros(size(b),'gpuArray');
    t_old = 1;
end

old_val = 0;
    
while(1)
   
    %% Step 1
    sol = Pc(b - lambda.*L(r,s));
    [tmpr,tmps] =  Lt(sol);
    
    r = r + (1./(8* lambda)).*tmpr;
    s = s + (1./(8* lambda)).*tmps;
    
    %Projection step Pp
    weights = max(1, sqrt(r.^2+s.^2));
    p = r./weights;
    q = s./weights;
    
    %% Step 2
    t = (1+sqrt(1+4*t_old^2))/2;
    t_old = t;
        
    %% Step 3
    r = p + ((t_old-1)/t) * (p-p_old);
    s = q + ((t_old-1)/t) * (q-q_old);
    
    p_old = p;
    q_old = q;
    
    %% Stopping crit
    tmp = (b-sol).^2;
    cur_val = .5*nansum(tmp(:)) +   lambda .* nansum(sqrt(tmpr(:).^2 + tmps(:).^2));
    rel_obj = abs(cur_val-old_val)/cur_val;

    if rel_obj<params.tvtol
        break;
    end

    old_val = cur_val;
end

if (params.usegpu==0)
    solution = sol;
else
    solution = gather(sol);
end

end


function [X] = fista_backtrack(grad, proj, Xinit, step, alpha, opts, calc_H, calc_G, constants)  
% function [X, iter, min_cost] = fista_general(grad,proj, Xinit, L, opts, calc_F)   
% Implementation of relaxed FISTA. 
% Ma, Yanting, et al. "Accelerated Image Reconstruction for Nonlinear Diffractive Imaging." 
% arXiv preprint arXiv:1708.01663 (2017).
%
% Solve the problem: x = arg min_x f(x) + lambda*h(x) 
%  Inputs:
%     grad   : a function calculating gradient of f(X) given X.
%     proj   : a function calculating pL(x) -- projection
%     Xinit  :initial guess.
%     step: Initial step size
%     alpha: FISTA relaxation
%       opts   : a struct
%           opts.lambda  : a regularization parameter, can be either a scalar or
%                           a weighted matrix.
%           opts.max_iteration: maximum iterations of the algorithm. 
%                           Default 300.
%           opts.tol     : a tolerance, the algorithm will stop if difference 
%                           between two successive X is smaller than this value. 
%                           Default 1e-8.
%           opts.verbose : showing f(x) after each iteration or not. 
%                           Default false. 
%       calc_G: a function calculating value of f(x) (data fidelity)
%       calc_H: a function calculating h(x) (regularizer)

    if ~isfield(opts, 'max_iteration')
        opts.max_iteration = 1100;
    end       
    if ~isfield(opts, 'tol')
        opts.tol = 2.4e-6;
    end
    if ~isfield(opts, 'verbose')
        opts.verbose = false;
    end

    x_old = Xinit;
    y_old = Xinit;
    t_old = 1;
    iter = 0;
    cost(1) = feval(calc_G,Xinit,constants) + feval(calc_H,Xinit,constants);
    lambda = opts.lambda;
    
    opts_proj = opts;
    opts_proj.lambda = lambda*step;

    while  iter < opts.max_iteration
        iter = iter + 1;

        yoldgrad = feval(grad, y_old, constants);
        x_new = feval(proj, y_old - (step)*yoldgrad, opts_proj);
        t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
        y_new = x_new + alpha*(t_old - 1)/t_new * (x_new - x_old);
        
        c = feval(calc_G, x_new,constants);
        oldc = feval(calc_G, y_old, constants);
        
        %% Line search
        while c > oldc + yoldgrad(:)'*(x_new(:)-y_old(:))+(1/(step*2))*norm((x_new-y_old))^2
            step = step*0.3;
            opts_proj.lambda = lambda*step;
            
            x_new = feval(proj, y_old - (step)*yoldgrad, opts_proj);
%             t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
            y_new = x_new + alpha*(t_old - 1)/t_new * (x_new - x_old);
            
            c = feval(calc_G, x_new, constants);
        end
        
        % Check the stop criteria
        e = norm(x_new - x_old)/numel(x_new);
        if iter>10
        if e < opts.tol
            break;
        end
        end
        
        %% Update
        x_old = x_new;
        t_old = t_new;
        y_old = y_new;
        
        cost(iter+1) = c + feval(calc_H, x_new,constants);
        
%         %% show progress
%         if opts.verbose
%         if ~mod(iter,50)
%                 disp(['iter = ',num2str(iter),' step = ',num2str(step),'  cost = ',num2str(cost(iter+1)),'  e = ',num2str(norm(x_new - x_old)/numel(x_new))]);
%         end 
%         end
        
        step = step*5; %Increasing step size slows each iterations but increases overall convergence speed
        
    end
    X = x_new;
end 
