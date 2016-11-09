module Examples

function generalized_rosenbrock(n, T = Float64)
  if n % 2 != 0
    error("only even dimensions are allowed")
  end
  F(x) = begin
    Y = Array{T}(n)
    for i in 1 : 2 : n
      Y[i] = T(1) - x[i]
    end
    for i in 2 : 2 : n
      Y[i] = T(10) * (x[i] - x[i-1]*x[i-1])
    end
    Y
  end
  J(x) = begin
    Y = zeros(T, n, n)
    for i in 1 : 2 : n
      Y[i, i] = T(-1)
    end
    for i in 2 : 2 : n
      Y[i, i] = T(10)
    end
    for i in 2 : 2 : n
      Y[i, i-1] = T(-20) * x[i-1]
    end
    Y
  end
  X0 = zeros(T, n)
  F , J, X0
end

function set_1(T = Float64)
  F(x) = [x[1] * x[1] + x[2] * x[2] - T(25),
          x[1] * x[1] - x[2] - T(2)]
  
  J(x) = [T(2) * x[1]  T(2) * x[2];
          T(2) * x[1]  T(-1)       ]
  
  X0 = [ T(1), T(1) ]
  
  F, J, X0
end

function set_2(T = Float64)
  F(x) = [x[1] + x[2] + x[3],
          x[1] * x[1] + x[2] * x[2] + x[3] *x[3] - T(2),
          x[1] * (x[2] + x[3]) + T(1)
         ]
  J(x) = [1           1         1;
          T(2)*x[1]   T(2)*x[2] T(2)*x[3];
          x[2] + x[3] x[1]      x[1]
         ]

  X = [T(3/4), T(1/2), T(-1/2)]

  F, J, X
end

function set_3(T = Float64)
  f1(x,y,z) = x*y-z^2-T(1)
  f2(x,y,z) = x*y*z+y^2-x^2-T(2)
  f3(x,y,z) = exp(x)+z-exp(y)-T(3)

  F(x) = [ f1(x[1],x[2],x[3]); f2(x[1],x[2],x[3]); f2(x[1],x[2],x[3]) ]

  J(x) = [ x[2]             x[1]              T(-2)x[3] ;
          x[2]*x[3]-T(2)x[1]  x[1]*x[3]+T(2)x[2]   x[1]*x[2] ;
          exp(x[1])        -exp(x[2])        T(1)
         ]
  F, J, ones(T, 3)
end

function set_4(T = Float64)
  F(x) = [x[1]*x[1], x[2]*x[2]]
  J(x) = [  T(2)*x[1] T(0)
          ; T(0)      T(2)*x[2]
         ]
  F, J, [T(1), T(1)]
end

function set_5(T = Float64)
  F(x) = [x[1]*x[1]*x[1] + T(7)*x[2] - T(7), x[1]*x[1] + x[2]*x[2] - T(1)]
  J(x) = [  T(3)*x[1]*x[1] T(7)
          ; T(2)*x[1]      T(2)*x[2]
         ]
  X0 = [T(1), T(2)]
  F, J, X0
end

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

function vandermonde(n)
  F(X) = begin
    B = Array{eltype(X)}(n)
    for i in 1 : n
      sum = 0
      for j in 1 : n
        sum += X[j] * (i + 1) ^ (j - 1)
      end
      B[i] = sum - ((i + 1) ^ n - 1) / i
    end
    B
  end
  J(X) = begin
    B = Array{eltype(X)}(n,n)
    for i in 1 : n
      for j in 1 : n
        B[i, j] = (i + 1) ^ (j - 1)
      end
    end
    B
  end
  (F, J)
end

function comparison_1(T = Float64)
  A = [ T(3.0) T(-13.0) T(9.0)   T(3.0)
      ; T(-6.0)   T(4.0) T(1.0) T(-18.0)
      ; T(6.0)  T(-2.0) T(2.0)   T(4.0)
      ; T(12.0)  T(-8.0) T(6.0)  T(10.0)
      ]
  B = [ T(-19), T(-34), T(16), T(26)]
  X₁ = gauss!(copy(A), copy(B))
  X₂ = gauss_naive!(copy(A), copy(B))
  X = [T(3), T(1), T(-2), T(1)]
  @printf("error for gauss_naive! : %.8e\nerror for gauss! : %.8e", 
          norm(X - X₂), norm(X - X₁))

end

function naive_fail(T = Float64, ϵ = T(1e-2))
  A = [ ϵ T(1)
      ; T(1)     T(1)
      ]
  
  B = [ T(1), T(2) ]

  A, B
end

function polynomial(n, T = Float64)
  p = 10 .* randn(n, n)
  b = zeros(T, n)
  F(x) = begin
    y = zeros(T, n)
    for i in 1 : n
      for j in 1 : n
        y[i] = y[i] + p[i, j] * x[j] * x[i]
        b[i] += p[i,j]
      end
      y[i] -= T(1000)
    end
    y
  end
  J(x) = begin
    y = zeros(T, n, n)
    for i in 1 : n
      for j in 1 : n
        y[i,j] += p[i,j]*x[i]
      end
      for k in 1 : n
        y[i,i] += p[i,k]*x[k]
      end
    end
    y
  end
  F, J, b
end
end # Examples
