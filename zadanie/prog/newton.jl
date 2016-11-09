include("utils.jl")
include("gauss.jl")

function newton(F, J, X
                ; solve = gauss!
                , ϵs = 1e-15
                , maxiter = 20
                , log = false
                , latex = false)
  ϵ = 1.0
  i = 0
  if log
    @printf("Iteration\t||F(X)||\t\tϵ\n")
  elseif latex
    @printf "Iteracja & \$\\|F(X)\\|\$ & \$\\epsilon\$ \\\\ \\hline\n"
  end
  while ϵ > ϵs && i < maxiter
    JX = J(X)
    FX = F(X)
    δX = solve(JX, FX)
    ϵ  = norm(δX) / norm(X)
    if log
      @printf("%d\t\t%e\t%e\n", i, norm(FX), ϵ)
    elseif latex
      @printf "%d & %e & %e \\\\ \\hline\n" i norm(FX) ϵ
    end
    X  = X - δX
    i += 1
  end
  if log
    @printf("\nFinal approximation:\n")
    print_array(X)
  end
end
