function print_array(A)
  @printf("[ ")
  for a in A
    @printf("%.3e ", a)
  end
  @printf("]\n")
end
