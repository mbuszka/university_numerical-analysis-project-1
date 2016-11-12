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
  elseif n == 6
    F, J, X0 = Examples.set_6()
  end
  @printf "test set %d begin\n" n
  X = newton(F, J, X0,
             solve   = gauss!,
             verbose = !latex,
             latex   = latex,
             ϵs      = 1e-12)
  @printf "end\n\n\n"
end

function run_set(n, log_friendly = false)
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
  elseif n == 6
    F, J, X0 = Examples.set_6()
  end
  Xs, iters, exs, efs = newton(F, J, X0, log = true, verbose = true)
  if log_friendly
    for i in 1:iters
      if exs[i] == 0
        exs[i] = eps(Float64)
      end
      if efs[i] == 0
        efs[i] = eps(Float64)
      end
    end
  end
  iters, exs[1:iters], efs[1:iters]
end

function dispatch_example(str, size, δ, T = Float64)
  if     str == "sin_cos"
    F, J, X0 = Examples.sin_cos(size, δ, T)
  elseif str == "sin_cos_2"
    F, J, X0 = Examples.sin_cos_2(size, δ, T)
  elseif str == "polynomial"
    F, J, X0 = Examples.polynomial(size, δ, T)
  elseif str == "vandermonde"
    F, J, X0 = Examples.vandermonde(size, δ, T)
  end
  F, J, X0
end

function calculate_accuracy(T)
  if     T == Float32
    return Float32(1e-8)
  elseif T == Float64
    return 1e-12
  elseif T == BigFloat # assume BigFloat precision 256 bit
    return big"1e-30"
  end
end

function show_generated(str; size = 8, T = Float64, latex = false)
  @printf "test generated function %s of size %d begin\n" str size
  F, J, X0 = dispatch_example(str, size, T)
  X, iter, exs, efs = newton(F, J, X0, 
                           solve   = gauss!,
                           log     = true,
                           verbose = !latex,
                           ϵs      = calculate_accuracy(T),
                           latex   = latex, 
                           maxiter = 5size)
  @printf "end\n\n\n"
end

function convergence_analysis(str, δ; 
                              verbose = false, 
                              turns   = 1000,
                              maxsize = 20,
                              T = Float64) 
  iters, divgs = Array{T}(maxsize-1), Array{T}(maxsize-1)
  acc = calculate_accuracy(T)
  if verbose
    @printf "convergence analysis of %s\n, δ = %.3f" str δ
    @printf "size\tavg iterations\tdivergent\n"
  end
  for i in 2 : maxsize
    iterm, divg = 0, 0
    for j in 1 : turns
      F, J, X0 = dispatch_example(str, i, δ, T)
      X, iter, exs, efs = newton(F, J, X0,
                                 solve   = gauss!,
                                 log     = true,
                                 ϵs      = acc,
                                 maxiter = 5i)
      if exs[iter] < acc
        iterm += iter
      else
        divg += 1
      end
    end
    iters[i-1] = float(iterm) / float(turns - divg)
    divgs[i-1] = float(divg)/float(turns)
    if verbose
      @printf("%d\t%.3f\t\t%.3f\n", i, iters[i-1], divgs[i-1])
    end
  end
  iters, divgs
end



