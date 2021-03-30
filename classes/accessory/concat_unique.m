function ab = concat_unique(a,b,epsilon)
    % Combine two arrays by comparing values and removing matching values
    % A matching value is who's difference from any other is less than some
    % defined epsilon
    % For matching values, use their mean.
    
    if nargin < 3
        epsilon = 1e-3;
    end
    
    ab = [];
    idx_ignore_b = [];
    
    % For each value in A, check if there are matching values in B. For all
    % matching values, append to AB the mean of the matching values. Then
    % add the indices of the matching values in B to a list to remember to
    % not consider them later.
    for ja = 1: length(a)
        idxb_matching = find( abs( b- a(ja)) < epsilon);
        idx_ignore_b = [idx_ignore_b, idxb_matching];
        ab = [ab, mean( [a(ja), b(idxb_matching)])];
    end
    
    if isempty(idx_ignore_b)
        ab = [ab,b];
        return
    end
    
    for jb = 1: length(b)
        if ~any(idx_ignore_b == jb)
            ab = [ab, b(jb)];
        end
    end
    
end