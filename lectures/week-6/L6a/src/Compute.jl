
# -- PRIVATE METHODS BELOW HERE ------------------------------------------------------------------------------- #
function _flux(problem::MyOptimalOpenExtentProblemCalculationModel)

    # initialize -
    results = Dict{String,Any}()
    c = problem.objective; # objective function coefficients

    # bounds -
    fluxbounds = problem.fluxbounds;
    lb = fluxbounds[:, 1]; # lower bounds
    ub = fluxbounds[:, 2]; # upper bounds
    
    # species bounds -
    speciesbounds = problem.speciesbounds;
    b = speciesbounds[:, 1]; # right-hand side vector
    
    # species constraints -
    A = problem.S; # constraint matrix

    # how many variables do we have?
    d = length(c);

    # Setup the problem -
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i,1] <= x[i=1:d] <= ub[i,1], start=0.0) # we have d variables
    
    # set objective function -   
    @objective(model, Max, transpose(c)*x);
    @constraints(model, 
        begin
            A*x >= -b # my material balance constraints 
        end
    );

    # run the optimization -
    optimize!(model)

    # check: was the optimization successful?
    @assert is_solved_and_feasible(model)

    # populate -
    x_opt = value.(x);
    results["argmax"] = x_opt
    results["objective_value"] = objective_value(model);

    # return -
    return results
end

function _flux(problem::MyPrimalFluxBalanceAnalysisCalculationModel)

    # initialize -
    results = Dict{String,Any}()
    c = problem.objective; # objective function coefficients

    # bounds -
    fluxbounds = problem.fluxbounds;
    lb = fluxbounds[:, 1]; # lower bounds
    ub = fluxbounds[:, 2]; # upper bounds
        
    # species constraints -
    A = problem.S; # constraint matrix

    # how many variables do we have?
    d = length(c);

    # Setup the problem -
    model = Model(GLPK.Optimizer)
    @variable(model, lb[i,1] <= x[i=1:d] <= ub[i,1], start=0.0) # we have d variables
    
    # set objective function -   
    @objective(model, Max, transpose(c)*x);
    @constraints(model, 
        begin
            A*x == 0 # my material balance constraints 
        end
    );

    # run the optimization -
    optimize!(model)

    # check: was the optimization successful?
    @assert is_solved_and_feasible(model)

    # populate -
    x_opt = value.(x);
    results["argmax"] = x_opt
    results["objective_value"] = objective_value(model);

    # return -
    return results
end
# -- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------- #

# -- PUBLIC METHODS BELOW HERE -------------------------------------------------------------------------------- #
"""
    softmax(x::Array{Float64,1})::Array{Float64,1}

Compute the softmax of a vector. 
This function takes a vector of real numbers and returns a vector of the same size with the softmax of the input vector, 
i.e., the exponential of the input vector divided by the sum of the exponential of the input vector.


### Arguments
- `x::Array{Float64,1}`: a vector of real numbers.

### Returns
- `::Array{Float64,1}`: a vector of the same size as the input vector with the softmax of the input vector.
"""
function softmax(x::Array{Float64,1})::Array{Float64,1}
    
    # compute the exponential of the vector
    y = exp.(x);
    
    # compute the sum of the exponential
    s = sum(y);
    
    # compute the softmax
    return y / s;
end

"""
    binary(S::Array{Float64,2})::Array{Int64,2}
"""
function binary(S::Array{Float64,2})::Array{Int64,2}
    
    # initialize -
    number_of_rows = size(S, 1);
    number_of_columns = size(S, 2);
    B = zeros(Int64, number_of_rows, number_of_columns);

    # main -
    for i ∈ 1:number_of_rows
        for j ∈ 1:number_of_columns
            if (S[i, j] != 0.0) # if the value is not zero, then B[i,j] = 1
                B[i, j] = 1;
            end
        end
    end

    # return -
    return B;
end

"""
    solve(model::AbstractFluxCalculationModel)
"""
function solve(model::AbstractFluxCalculationModel)
    return _flux(model);
end

# -- PUBLIC METHODS ABOVE HERE -------------------------------------------------------------------------------- #