function [ out ] = iidbernoulli( sz, theta )
%IIDBERNOULLI   iid Bernoulli sparsity pattern w/ at least 1 nonzero.
    atleast1sparse = false;
    while ~atleast1sparse
        out = double(rand(sz) <= theta);
        atleast1sparse = sum(out(:) ~= 0) > 0;
    end
end

