function [ A, X, info ] = AXsolve_PALM( Y, k, lambda, varargin )
%AXSOLVE_PALM   Solve for A and X in the SBD objective with PALM.
    addpath([fileparts(mfilename('fullpath')) '\helpers']);
    m = size(Y);    
    m = m(1:2); 
    n = size(Y,3);
    
    % Default parameter values
    d_maxit = 1e3;
    d_Leps = 1.1;
    d_gradtol = 1e-3;
    d_verbose = floor(d_maxit / 10);
    d_showvars = false;
    d_dispfun = @() 0;
    
    p = inputParser;    
    addParameter(p, 'A0', proj2oblique(randn([k n])), @(A) (ismatrix(A) || ndims(A)==3));
    addParameter(p, 'X0', randn(m), @(X) ismatrix(X));
    
    addParameter(p, 'maxit', d_maxit, @(x) x>=0);
    addParameter(p, 'Leps', d_Leps, @(x) x>=0);
    addParameter(p, 'gradtol', d_gradtol, @(x) x>=0);
    addParameter(p, 'verbose', d_verbose, @isnumeric);
    addParameter(p, 'showvars', d_showvars);
    addParameter(p, 'dispfun', d_dispfun);
    
    parse(p, varargin{:});
    A = p.Results.A0;
    X = p.Results.X0;
    p = rmfield(p.Results, {'A0', 'X0'});
    
    %% Iterate:    
    R = zeros([m n]);       % the residual (\mcA \boxast X - \mcY)
    GX = zeros([m n]);      % gradient information for X
    GA = zeros([k n]);      % gradient information for A
    Ahat = zeros([m n]);    % FFT of A;
    
    repeat = true;
    zeroX = false;          % watch for zero-solution
    costs = NaN([p.maxit 2]);
    it = 0;
    while repeat
        it = it+1;
        
        % Update X
        for i = 1:n
            R(:,:,i) = convfft2(A(:,:,i), X) - Y(:,:,i);
            GX(:,:,i) = convfft2(A(:,:,i), R(:,:,i), true);
            Ahat(:,:,i) = fft2(A(:,:,i), m(1), m(2));
        end
        costs(it, 1) = norm(R(:)-Y(:))^2/2 + lambda*norm(X(:),1);
        
        tmp = sum(abs(Ahat).^2, 3);
        LX = p.Leps * max(tmp(:));
        
        gradX = sum(GX, 3);
        X = X - gradX/LX;                       % gradient step
        X = sign(X).*max(abs(X)-lambda/LX, 0);  % prox. step (soft-threshold)
        
        if ~norm(X(:),1)
            zeroX = true;
            warning('X == zeros: output is invalid. Try lowering lambda.');
        end
            
        % Update A
        for i = 1:n
            R(:,:,i) = convfft2(A(:,:,i), X) - Y(:,:,i);            
            tmp = convfft2(X, R(:,:,i), true, m+k-1, m);
            GA(:,:,i) = tmp(1:k(1), 1:k(2));
        end
        costs(it, 2) = norm(R(:)-Y(:))^2/2 + lambda*norm(X(:),1);
        
        tmp = abs(fft2(X));
        LA = p.Leps * max(tmp(:))^2;
        
        A = A - GA/LA;          % gradient step;
        A = proj2oblique(A);    % prox. step (projection)
        
        % Check for repeat conditions
        gnormX = norm(gradX,'fro')/sqrt(numel(X));
        gnormA = norm(GA(:))/sqrt(numel(A));
        gtolreached = prod([gnormX gnormA] <= p.gradtol);
        repeat = (it < p.maxit) && ~zeroX && ~gtolreached;
    end
    
    info.it = it;
    info.maxitreached = it >= p.maxit;
    info.zeroX = zeroX;
    info.gradtolreached = gtolreached;
    info.vars = {'X' 'A'};
    info.gradnorm = [gnormX gnormA];
    info.L = [LX LA];
    info.costs = costs;     % note that these are costs BEFORE update
end


%{
    %% Initialize arguments
    oargnum = 4;
    if nargin < oargnum-1
        error(['Arguments up to argument ' num2str(oargnum-1) ' are mandatory.']);
    end
        
    if (nargin < oargnum) || isempty(params)
        params = 1e3;
    end; oargnum = oargnum+1;
    if (nargin < oargnum) || isempty(A0)
        A = proj2oblique(randn([k n]));
    else
        A = A0;
    end; oargnum = oargnum+1;
    if (nargin < oargnum) || isempty(X0)
        X = randn([m n]);
    else
        X = X0;
    end
%}