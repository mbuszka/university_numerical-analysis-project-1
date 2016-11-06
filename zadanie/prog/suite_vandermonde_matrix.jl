function vandermonde_matrix(n, T = Float64)
  A = Array{T}(n,n)
  B = Array{T}(n)
  for i in 1 : n
    for j in 1 : n
      A[i, j] = (i + 1) ^ (j - 1)
    end
    B[i] = ((i +1) ^ n - 1) / i
  end
  A, B
end

function suite_vandermonde_matrix(T)
  range = if T == Float32
        5 : 10
      elseif T == Float64
        8 : 15
      elseif T == BigFloat
        15 : 25
      end
  for n in $range
    A, B = vandermonde_matrix(n, $T)
    X    = gauss!(A, B)
    ϵ    = norm(X - ones($T, n))
    @printf("Error for matrix size %d is %e", n, ϵ)
  end
end

function vandermonde_test(solver, nmax, T = Float64)
  for n in 2 : nmax
    A, B = vandermonde_matrix(n, T)
    X = solver(A, B)
    ϵ = norm(X - ones(eltype(X), n))
    @printf("Error for matrix size n = %d is ϵ = %.8e\n", n, ϵ)
  end
end

function comparison1(T = Float64)
  X₁ = gauss!([ T(3.0) T(-13.0) T(9.0)   T(3.0)
              ;T(-6.0)   T(4.0) T(1.0) T(-18.0)
              ; T(6.0)  T(-2.0) T(2.0)   T(4.0)
              ;T(12.0)  T(-8.0) T(6.0)  T(10.0)],
              [ T(-19), T(-34), T(16), T(26)])
  X₂ = gauss_naive!([ T(3.0) T(-13.0) T(9.0)   T(3.0)
              ;T(-6.0)   T(4.0) T(1.0) T(-18.0)
              ; T(6.0)  T(-2.0) T(2.0)   T(4.0)
              ;T(12.0)  T(-8.0) T(6.0)  T(10.0)],
              [ T(-19), T(-34), T(16), T(26)])
  X = [T(3),T(1), T(-2), T(1)]
  @printf("error for gauss_naive! : %.8e\nerror for gauss! : %.8e", 
          norm(X - X₂), norm(X - X₁))

end
