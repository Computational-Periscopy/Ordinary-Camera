function [reconstruction] = reconstruct_tv(Ar,Ag,Ab, measurement, reg, NumBlocks)
%% Forms the reconstruction given a measurement image and using the A matrices.
% Includes non-negativity constraint

    
    mr = measurement(:,:,1);
    mg = measurement(:,:,2);
    mb = measurement(:,:,3);    
    
    % Initial estimates
    tik_reg_param =  0; %100000;
    testr = (Ar'*Ar + tik_reg_param*size(Ar,1)*(ones(size(Ar,2))))\(Ar') * mr(:);
    testg = (Ag'*Ag + tik_reg_param*size(Ag,1)*(ones(size(Ag,2))))\(Ag') * mg(:);
    testb = (Ab'*Ab + tik_reg_param*size(Ab,1)*(ones(size(Ab,2))))\(Ab') * mb(:);
    
%     testr(testr<0) = 0;
%     testb(testb<0) = 0;
%     testg(testg<0) = 0;
%   
%     testr(testr>1) = 1;
%     testb(testb>1) = 1;
%     testg(testg>1) = 1;
    
    %%
    ftv.params.constraint = 'non-neg'; %Positive values only

    testr = co_prox_tv(reshape(testr(:),NumBlocks(1),NumBlocks(2)), reg(1), ftv.params);
    testg = co_prox_tv(reshape(testg(:),NumBlocks(1),NumBlocks(2)), reg(2), ftv.params);
    testb = co_prox_tv(reshape(testb(:),NumBlocks(1),NumBlocks(2)), reg(3), ftv.params);
     
    reconstruction(:,:,1) = testr;
    reconstruction(:,:,2) = testg;
    reconstruction(:,:,3) = testb;


end

function [solution] = co_prox_tv(input, lambda, params)
%% Solves prox_tv(Y) = argmin_X  1/2 ||X-Y||^2 + lambda*TV(X)
% Using method from Beck & Teboulle (2009)
% TV = Isotropic TV with reflexive boundary conditions, i.e d/dx(x>=size(input)) = 0;

% input = input matrix/image Y
% Parameters
%    params.iterations  => Perform number of iterations 
%    params.tol         => Stopping tolerance, stops when:
%                          (error(i)-error(i-1))/error(i)<param.tol
%                          Defaults to 1e-6
%    params.constraint  => 'none' - default
%                          'non-neg'
%    params.usegpu      => Use GPU arrays. 0 = No, 1 = Yes. Default = 0;

%% Optional inputs
if nargin==2
    params.tol = 1e-6;
    params.constraint = 'none';
    params.usegpu = 0;
else
    if ~isfield(params, 'tol') && ~ isfield(params, 'iterations')
        params.tol = 1e-6;
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
    Pc = @(x) x.*(x>0);  %Project onto constrained set i.e non-negatvity contraints 
%     Pc = @(x) min(1,max(0,x));
elseif strcmp(params.constraint,'bounds')
    Pc = @(x) max(min(1,x),0);  %Project onto constrained set i.e non-negatvity contraints 
%     Pc = @(x) min(1,max(0,x));
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
    
if isfield(params, 'tol')  %Solve using stopping tolerance
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

    if rel_obj<params.tol
        break;
    end

    old_val = cur_val;
end
else
for it=1:params.iterations
   
    %% Step 1
    sol = Pc(b -  lambda.*L(r,s));
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
    
end  %Solve using set number of iterations
end

if (params.usegpu==0)
    solution = sol;
else
    solution = gather(sol);
end


end


