function [problem] = fnInitMethod(sMethod)

problem = struct();
problem.options.method = sMethod;

switch problem.options.method
    case 'trapezoid'
        problem.options.trapezoid.nGrid = 70;
    case 'hermiteSimpson'
        problem.options.hermiteSimpson.nSegment = 155;
    otherwise
        error('Invalid method.');
end

end