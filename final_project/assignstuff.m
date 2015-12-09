function varargout = assignstuff(n)
varargout = cell(1,length(n));
for i = 1:length(n)
    varargout{i} = n(i);
end