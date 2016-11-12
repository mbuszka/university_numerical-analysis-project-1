include("utils.jl")
include("gauss.jl")

function newton(F, J, X
                ; solve::Function = gauss!
                , ϵs = 1e-12
                , maxiter::Int64 = 20
                , log::Bool = false
                , verbose::Bool = false
                , latex::Bool = false)
  T = eltype(X)
  ϵf, ϵx = T(1), T(1)
  i = 0
  if log
    exs = zeros(T, maxiter)
    efs = zeros(T, maxiter)
    Xs  = zeros(T, size(X, 1), maxiter)
  end
  if verbose
    @printf "Iteration\t||F(X)||\t\tϵ\n"
  elseif latex
    @printf "Iteracja & \$\\|F(X)\\|\$ & \$\\epsilon\$ \\\\ \\hline\n"
  end
  JX = J(X)
  FX = F(X)
  while (ϵx > ϵs || ϵf > ϵs) && i < maxiter
    i += 1
    δX = solve(JX, FX)
    ϵx = norm(δX) / norm(X)
    X  = X - δX
    FX = F(X)
    JX = J(X)
    ϵf = norm(FX)
    if log
      Xs[:, 1] = X
      exs[i] = ϵx
      efs[i] = ϵf
    end
    if verbose
      @printf "%d\t\t%e\t%e\n" i ϵf ϵx
    elseif latex
      @printf "%d & %e & %e \\\\ \n" i ϵf ϵx
    end
  end
  if verbose
    @printf "\nFinal approximation:\n"
    print_array(X)
  end
  if log
    return Xs, i, exs, efs
  else
    return X
  end
end
