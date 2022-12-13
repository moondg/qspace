function i=isfield(A,f)
% wrapper for builtin isfield
% Wb,Feb26,16

  i=isfield(struct(A),f);

end

