function [problem] = fnInitMethod(sMethod,nOfPoints,...
    optimalDesignFlag,fuelSavingFlag,tireEnergySavingFlag)

problem = struct();
problem.options.method = sMethod;

switch problem.options.method
    case 'trapezoid'
        problem.options.trapezoid.nGrid = nOfPoints;
    case 'hermiteSimpson'
        problem.options.hermiteSimpson.nSegment = 350;
    otherwise
        error('Invalid method.');
end

% optimalDesignFlag = false;
% fuelSavingFlag    = true;

problem.options.optimalDesignFlag       = optimalDesignFlag; 
problem.options.fuelSavingFlag          = fuelSavingFlag; 
problem.options.tireEnergySavingFlag    = tireEnergySavingFlag; 

end