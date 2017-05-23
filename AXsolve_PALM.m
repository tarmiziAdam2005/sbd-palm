function [ A, X, info ] = AXsolve_PALM( Y, k, lambda, params, A0, X0 )
%AXSOLVE_PALM   Solve for A and X in the SBD objective with PALM.
    addpath([fileparts(mfilename('fullpath')) '\helpers']);

    m = size(Y);    
    m = m(1:2); 
    n = size(Y,3);

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
    
    %% Parameters
    maxit = params;        % maximum iterations
    Leps = 1.1;         % Lipschitz safety factor (>1).
    
    %% Iterate:    
    R = zeros([m n]);       % the residual (\mcA \boxast X - \mcY)
    GX = zeros([m n]);      % gradient information for X
    GA = zeros([k n]);      % gradient information for A
    Ahat = zeros([m n]);    % FFT of A;
    
    repeat = true;
    zeroX = false;          % watch for zero-solution
    costs = NaN([maxit 2]);
    it = 1;
    while repeat
        % Update X
        for i = 1:n
            R(:,:,i) = convfft2(A(:,:,i), X) - Y(:,:,i);
            GX(:,:,i) = convfft2(A(:,:,i), R(:,:,i), true);
            Ahat(:,:,i) = fft2(A(:,:,i), m(1), m(2));
        end
        costs(it, 1) = norm(R-Y,'fro')^2/2 + lambda*norm(X(:),1);
        
        tmp = sum(abs(Ahat).^2, 3);
        LX = Leps * max(tmp(:));
        
        X = X - sum(GX, 3)/LX;                  % gradient step
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
        costs(it, 2) = norm(R-Y,'fro')^2/2 + lambda*norm(X(:),1);
        
        tmp = abs(fft2(X));
        LA = Leps * max(tmp(:))^2;
        
        A = A - GA/LA;          % gradient step;
        A = proj2oblique(A);    % prox. step (projection)
        
        % Check for repeat conditions
        repeat = (it <= maxit) && ~zeroX;
        it = it+1;
    end
    
    info = costs;
end


