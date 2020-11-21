function fitOpt = getFitOptions(solver, init)

if ~isempty(solver)
    fitOpt.solver = solver;
end

if ~isempty(init)
    fitOpt.init = init;
end

end

