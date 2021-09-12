function out = boolNaN(bool)
    out = zeros(size(bool));
    bool = logical(bool);
    out(bool) = NaN(nnz(bool),1);
end