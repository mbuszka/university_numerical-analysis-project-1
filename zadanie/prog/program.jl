function newton(F, J, X, solve = (A, Y) -> \(A, Y), ϵs = 1e-15, maxiter = 20, log = false)
  ϵ = 1.0
  i = 0
  while ϵ > ϵs && i < maxiter
    JX = J(X)
    FX = F(X)
    δX = solve(JX, FX)
    # @show JX
    # @show δX
    ϵ  = norm(δX) / norm(X)
    # @show X
    if log
      @printf("ϵ = %e\n", ϵ)
    end
    X  = X - δX
    i += 1
  end
  @printf("final approximation :\n  ")
  @show X
  @show F(X)
end


"""
Oblicza rozwiązanie układu równań AX = B metodą eliminacji Gaussa

Korzysta z naiwnego algorytmu zaadaptowanego z książki Numerical Mathematics
and Computing Warda Cheney'a i Davida Kincaida
"""
function gauss_naive!(A, Y)
  n = length(Y)
  X = Array{eltype(A)}(n)
  for k in 1 : n-1
    for i = k + 1 : n
      mult = A[i, k] / A[k, k]
      A[i, k] = 0 
      for j in k + 1 : n
        A[i, j] -= mult * A[k, j]
      end
      Y[i] -= mult * Y[k]
    end
  end
  X[n] = Y[n] / A[n, n]
  for i in n-1 : -1 : 1
    sum = Y[i]
    for j in i + 1 : n
      sum -= A[i, j] * X[j]
    end
    X[i] = sum / A[i, i]
  end
  X
end


"""
Oblicza rozwiązanie układu równań metodą eliminacji Gaussa.

Korzysta z algorytmu skalowanego elementu głównego, przekształca
macierz wejściową A, tak że w górnej trójkątnej jej części znajduje się
jej postać schodkowa, natomiast w dolnej zapisane są użyte mnożniki
które wykorzystywane są w kroku podstawiania

Algorytm został zaadaptowany z książki Numerical Mathematics and Computing
Warda Cheney'a i Davida Kincaida
"""
function gauss!(A, Y)
  n = length(Y)
  l = gauss_eliminate!(A, n)
  X = gauss_substitute!(A, Y, l, n)
  X
end


"""
Implementuje krok eliminacji, zapisując mnożniki w dolnej trójkątnej
części macierzy A i zwracając permutację wierszy macierzy A
"""
function gauss_eliminate!(A, n)
  s = Array{eltype(A)}(n) # Scale vector
  l = Array{Int}(n)       # Index permutation
  for i in 1 : n
    l[i] = i
    smax = 0
    for j in 1 : n
      smax = max(smax, abs(A[i, j]))
    end
    s[i] = smax
  end
  for k in 1 : n-1
    rmax = 0
    for i in k : n
      r = abs(a[l[i], k] / s[l[i]])
      if r > rmax
        rmax = r
        j = i
      end
    end
    swap(l[j], l[k])
    for i = k + 1 : n
      mult = A[l[i], k] / A[l[k], k]
      A[l[i], k] = 0 
      for j in k + 1 : n
        A[l[i], j] -= mult * A[l[k], j]
      end
    end
  end
  l
end

"""
Implementuje krok podstawiania
"""
function gauss_substitute!(A, Y, l, n)
  X = Array{eltype(A)}(n) # Solution vector
  X[n] = Y[l[n]] / A[l[n], n]
  for i in n-1 : -1 : 1
    sum = Y[l[i]]
    for j in i + 1 : n
      sum -= A[l[i], j] * X[j]
    end
    X[i] = sum / A[l[i], i]
  end
  X
end

function test_naive(size = 10)  
  for n in 4 : size
    A = Array{Float32}(n,n)
    Y = Array{Float32}(n)
    for i in 1 : n
      for j in 1 : n
        A[i, j] = (i + 1) ^ (j - 1)
      end
      Y[i] = ((i +1) ^ n - 1) / i
    end
    T = \(A, Y)
    X = gauss_naive!(A, Y)
    @show X
    # @show T
  end
end

function vandermonde(n)
  F(X) = begin
    Y = Array{eltype(X)}(n)
    for i in 1 : n
      sum = 0
      for j in 1 : n
        sum += X[j] * (i + 1) ^ (j - 1)
      end
      Y[i] = sum - ((i + 1) ^ n - 1) / i
    end
    Y
  end
  J(X) = begin
    Y = Array{eltype(X)}(n,n)
    for i in 1 : n
      for j in 1 : n
        Y[i, j] = (i + 1) ^ (j - 1)
      end
    end
    Y
  end
  (F, J)
end

function test_vandermonde(n) 
  @printf("Testing on vandermonde-like system of equtions\n")
  for i in 2 : n
    (F, J) = vandermonde(i)
    @printf("number of equations : %d\n", i)
    newton(F, J, rand(i))
  end
end

f1(x,y,z) = x*y-z^2-1
f2(x,y,z) = x*y*z+y^2-x^2-2
f3(x,y,z) = exp(x)+z-exp(y)-3

F(x) = [ f1(x[1],x[2],x[3]); f2(x[1],x[2],x[3]); f2(x[1],x[2],x[3]) ]

J(x) = [ x[2]             x[1]              -2x[3] ;
         x[2]*x[3]-2x[1]  x[1]*x[3]+2x[2]   x[1]*x[2] ;
         exp(x[1])        -exp(x[2])        1
       ]

newton(F, J, ones(Float64, 3))
newton(F, J, ones(Float64, 3), gauss_naive) 
test_vandermonde(15)
