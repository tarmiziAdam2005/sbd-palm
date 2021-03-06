function [ A, X, info ] = AXsolve_PALM( Y, k, lambda, varargin )
%AXSOLVE_PALM   Solve for A and X in the SBD objective with PALM.
    addpath([fileparts(mfilename('fullpath')) '\helpers']);
    m = size(Y);    
    m = m(1:2); 
    n = size(Y,3);
    
    % Default parameter values
    d_maxit = 1e3;
    d_Leps = 1.1;
    d_difftol = 1e-3;
    d_verbose = [];
    d_showvars = false;
    d_dispfun = @(it, Y, A, X, costs) 0;
    
    p = inputParser;    
    addParameter(p, 'A0', proj2oblique(randn([k n])), @(A) (ismatrix(A) || ndims(A)==3));
    addParameter(p, 'X0', randn(m), @(X) ismatrix(X));
    
    addParameter(p, 'maxit', d_maxit, @(x) x>=0);
    addParameter(p, 'Leps', d_Leps, @(x) x>=0);
    addParameter(p, 'difftol', d_difftol, @(x) x>=0);
    addParameter(p, 'verbose', d_verbose, @isnumeric);
    addParameter(p, 'showvars', d_showvars);
    addParameter(p, 'dispfun', d_dispfun);
    
    parse(p, varargin{:});
    A = p.Results.A0;
    X = p.Results.X0;
    p = rmfield(p.Results, {'A0', 'X0'});
    if isempty(p.verbose)
        p.verbose = floor(p.maxit/10);
    end
    
    %% Iterate:    
    R = zeros([m n]);       % the residual (\mcA \boxast X - \mcY)
    GX = zeros([m n]);      % gradient information for X
    rGA = zeros([k n]);     % Riemannian gradient information for A
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
        Xold = X;
        X = X - gradX/LX;                       % gradient step
        X = sign(X).*max(abs(X)-lambda/LX, 0);  % prox. step (soft-threshold)
        
        if ~norm(X(:),1)
            zeroX = true;
            warning('X == zeros: output is invalid. Try lowering lambda.');
        end
            
        % Update A
        for i = 1:n
            Ai = A(:,:,i);
            PA = Ai(:)*Ai(:)';
            R(:,:,i) = convfft2(Ai, X) - Y(:,:,i);            
            tmp2 = convfft2(X, R(:,:,i), true, m+k-1, m);
            GAi = tmp2(1:k(1), 1:k(2));
            rGA(:,:,i) = GAi - reshape(PA*GAi(:), size(GAi));
        end
        costs(it, 2) = norm(R(:)-Y(:))^2/2 + lambda*norm(X(:),1);
        
        tmp = abs(fft2(X));
        LA = p.Leps * max(tmp(:))^2;
        
        Aold = A;
        A = A - rGA/LA;         % Riemannian gradient step;
        A = proj2oblique(A);    % projection step;
        
        % Check for repeat conditions
        diffX = norm(X-Xold,'fro')/sqrt(numel(X));
        diffA = norm(A(:)-Aold(:))/sqrt(numel(A));
        gradnormA = norm(rGA(:))/sqrt(prod(k))/n;
        tolreached = prod([diffX diffA gradnormA] <= p.difftol);
        repeat = (it < p.maxit) && ~zeroX && ~tolreached;
        
        % Display update info:
        if p.verbose > 0 && (~mod(it, p.verbose) || it==1 || ~repeat)
            updatestr = sprintf('Iteration %s%d', ...
                blanks(floor(log10(p.maxit))-floor(log10(it))), it);
            if ~repeat
                updatestr = [updatestr ' (final)']; %#ok<*AGROW>
            end
            updatestr = [updatestr ...
                sprintf(': cost = %.4e, diff = %.4e, ||rgrad(A)||_F = %.4e.', ...
                costs(it,2), max(diffX, diffA), gradnormA)];
            disp(updatestr);
            
            p.dispfun(it, Y, A, X, costs(:,2));
        end
    end
    
    info.it = it;
    info.maxitreached = it >= p.maxit;
    info.zeroX = zeroX;
    info.tolreached = tolreached;
    info.tolcrit = {'diffX' 'diffA' 'gradnormA'};
    info.critvals = [diffX diffA gradnormA];
    info.L = [LX LA];
    info.costs = costs;     % note that these are costs BEFORE update
end