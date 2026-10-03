function [nInputs, inputNames] = getFuncHandleInputs(func_handle)
arguments
	func_handle
end
fh=func_handle;
% GETFUNCHANDLEINPUTS Return number of inputs and their names for a function handle
%
%   [nInputs, inputNames] = getFuncHandleInputs(fh)
%   returns the number of inputs (nInputs) and a string array of input
%   argument names (inputNames). Works for anonymous functions and file‑based
%   functions; built‑ins return inputNames = string.empty.

    s = functions(fh);

    if strcmp(s.type,'anonymous')
        header = regexp(s.function, '@\(([^)]*)\)', 'tokens','once');
        if isempty(header) || isempty(header{1})
            inputNames = string.empty;
        else
            inputNames = string(strtrim(strsplit(header{1}, ',')));
        end
        nInputs = numel(inputNames);
    else
        nInputs = nargin(fh);
        try
            fname = which(s.function);
            firstLine = regexp(fileread(fname), '^\s*function.*?\(([^)]*)\)', ...
                               'tokens','once','lineanchors');
            if isempty(firstLine) || isempty(firstLine{1})
                inputNames = string.empty;
            else
                inputNames = string(strtrim(strsplit(firstLine{1}, ',')));
            end
        catch
            inputNames = string.empty;
        end
	end
	inputNames=inputNames(:)';
end
