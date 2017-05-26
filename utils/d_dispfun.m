function d_dispfun( Y, A, X, costs, A0, X0, varargin )
%D_DISPFUN  Default display function for use when A0 and X0 are available.
%   sliceidx should be a function handle that generates natural numbers

    p = inputParser;
    addOptional(p, 'it', []);
    addOptional(p, 'sliceidx', @() 1);
    addOptional(p, 'abs', true);
    parse(p, varargin{:}); p = p.Results;
    i = p.sliceidx();

    subplot(231); imagesc(abs(Y(:,:,i))); 
    subplot(234); imagesc(abs(convfft2(A(:,:,i), X)));

    subplot(232); imagesc(abs(A0(:,:,i))); 
    subplot(235); imagesc(abs(A(:,:,i)));

    subplot(233); imagesc(abs(X0)); 
    subplot(236); imagesc(abs(X));

    drawnow;
end

