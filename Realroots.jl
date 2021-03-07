
using LinearAlgebra
using PolynomialRoots
 S=roots([2,-9,-1,8,5])
  A=Float64[0 ,0 ,0 ,0]
    for i=1:4 
     a=real(S[i])
     b=imag(S[i])
    if b != 0
        continue
        
    end
     A[i]=a
     #println(a)
     i+= 1
end
S= filter(!iszero, A)



