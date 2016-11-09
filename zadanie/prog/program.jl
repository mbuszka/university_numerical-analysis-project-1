include("newton.jl")
include("examples.jl")

function show_naive_fail(solve = gauss_naive!, latex = false)
  ϵ = 1e-1
  if latex
    @printf "\\epsilon & \$x_1\$ & \$x_2\$ \\\\ \\hline\n"
  else
    @printf "ϵ\t\tx₁\t\tx₂\n"
  end
  for i in 1 : 20
    A, B = Examples.naive_fail(Float64, ϵ)
    X = solve(A, B)
    if latex
      @printf "%.5e & %.7e & %.7e \\\\ \\hline\n" ϵ X[1] X[2]
    else
      @printf "%.5e\t%.7e\t%.7e\n" ϵ X[1] X[2]
    end
    ϵ = ϵ / 10
  end
end

function compare_gauss(latex = false)
  println("Naive Gauss implementation fails")
  show_naive_fail(gauss_naive!, latex)
  println("\nBut better one doesn't have a problem")
  show_naive_fail(gauss!, latex)
  println()
end

function show_set(n, latex = false)
  if     n == 1
    F, J, X0 = Examples.set_1()
  elseif n == 2
    F, J, X0 = Examples.set_2()
  elseif n == 3
    F, J, X0 = Examples.set_3()
  elseif n == 4
    F, J, X0 = Examples.set_4()
  elseif n == 5
    F, J, X0 = Examples.set_5()
  end
  X = newton(F, J, X0, solve = gauss!, log = !latex, latex = latex)
end
setprecision(256)
F, J, X = Examples.polynomial(3, BigFloat)
@show F([1.0,1.0,1.0])
@show F([1.0,1.0,1.0])
@show J([1.0,1.0,1.0])
newton(F, J, [BigFloat(1.2), BigFloat(0.8), BigFloat(0.9)], log = true, ϵs = big"1e-50")

compare_gauss()
show_set(5, false)
   
