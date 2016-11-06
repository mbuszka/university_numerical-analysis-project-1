include("gauss.jl")
include("newton.jl")
using Gauss, Newton

gauss_suites = ["vandermonde_matrix"]
newton_suites = ["3x3-1"]
setprecision(256)

function runnewton(T)
  for s in newton_suites
    include("suite_$s.jl")
    runsuite(T)
  end
end

function rungauss(T)
  for s in gauss_suites
    include("suite_$s.jl")

  end
end

for T in [Float32, Float64, BigFloat]
  rungauss(T)
  # runnewton(T)
end

macro runsuite(name, T)
  return :()
