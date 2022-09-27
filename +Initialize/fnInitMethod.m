function [problem] = fnInitMethod(sMethod,optimalDesignFlag,fuelSavingFlag)

problem = struct();
problem.options.method = sMethod;

switch problem.options.method
    case 'trapezoid'
        problem.options.trapezoid.nGrid = 500;
    case 'hermiteSimpson'
        problem.options.hermiteSimpson.nSegment = 155;
    otherwise
        error('Invalid method.');
end

% optimalDesignFlag = false;
% fuelSavingFlag    = true;

problem.options.optimalDesignFlag = optimalDesignFlag; 
problem.options.fuelSavingFlag    = fuelSavingFlag; 

end