using Roots
  using LinearAlgebra
  
  λ₁=2
  λ₂=1
  β=2
  h₁ =1
  h₂ =1
  f(p₁) = (((λ₁-λ₂)p₁+(λ₂*h₁))^2)*(λ₁*(p₁^2)-β)+((λ₁^2)*(λ₂)*(h₂^2)*(p₁^2)) 
   find_zeros(f, -10,  10)