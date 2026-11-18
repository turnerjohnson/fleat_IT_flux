function Ssize = whosstruct(S)
    for fld = fieldnames(S)'
        val = S.(fld{1});
        tmp = whos('val');
        Ssize.(fld{1}) = tmp.size;
    end
end