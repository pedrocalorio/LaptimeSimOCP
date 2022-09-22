function [problem] = fnInitFunctionHandles(problem)

problem.func.stateDynamics = @(x)( Controller.fnDynamics(x,u,problem.dsSystem.vd) );
problem.func.objective = @(x,u)( Controller.fnObjective(x,problem.dsSystem.td) );
problem.func.constraints = @(x)( Controller.fnConstraints(x,u) );

end