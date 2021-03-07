  using LinearAlgebra
  using PolynomialRoots
  #=λ₁=2
  λ₂=1
  β=1
  h₁ =3
  h₂ =0
  =#
#=============================================#
function Quadratic(x,A,b)
    return dot(x,A*x)+2*dot(b,x)
end
   

#====================================#
        uₖ =A*x₀+b
        vₖ =a-x₀
        
        uₖdotvₖ=dot(uₖ,vₖ)
        vₖdotvₖ=dot(vₖ,vₖ)
        normvₖ=norm(vₖ)
        normuₖ=norm(uₖ)
        zₖ =uₖ-((uₖdotvₖ)/(vₖdotvₖ))*vₖ
        normzₖ=norm(zₖ)
        zₖdotzₖ=dot(zₖ,zₖ)
        uₖdotzₖ=dot(uₖ,zₖ)
        Avₖ=A*vₖ 
        Azₖ=A*zₖ 
        vₖdotAvₖ=dot(vₖ,Avₖ)
        zₖdotAzₖ=dot(zₖ,Azₖ)
        vₖdotAzₖ=dot(vₖ,Azₖ)
        zₖdotAvₖ=dot(zₖ,Avₖ)

        Aₖ=[(vₖdotAvₖ)/ (normvₖ^2)  (vₖdotAzₖ)/(normvₖ*normzₖ) ; (zₖdotAvₖ)/(normvₖ*normzₖ)  (zₖdotAzₖ)/(normzₖ^2)]
        bₖ=[(uₖdotvₖ)/(normvₖ) ;(uₖdotzₖ)/(normzₖ)]

        D, Q = eigen(Aₖ)
       Qbₖ=Q*bₖ
       D⁻¹Qbₖ=D.\Qbₖ
       #Aₖ=Q'*D*Q
       β=dot(Qbₖ,D⁻¹Qbₖ)
       aₗ=[normvₖ,0]

       aₕ=Q*aₗ+D⁻¹Qbₖ
       #======================#
       h₁=aₕ[1]
       h₂=aₕ[2]
       λ₁=D[1]
       λ₂=D[2]

  #P=((λ₁*(λ₁-λ₂)^2))*(p₁)^4+(2*λ₁*λ₂*h₁*(λ₁-λ₂))*(p₁)^3+((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))*(p₁^2)-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)*(p₁)-(β*λ₂^2*h₁^2)
  a4=((λ₁*(λ₁-λ₂)^2))
  a3=(2*λ₁*λ₂*h₁*(λ₁-λ₂))
  a2=((λ₁*λ₂)*((λ₂*h₁^2)+(λ₁*h₂^2))-(β*(λ₁-λ₂)^2))
  a1=-2*β*λ₂*h₁*(λ₁-λ₂)
  a0=-1*(β*λ₂^2*h₁^2)
  S=roots([a0,a1, a2, a3, a4])
  #choose Real roots
  #=G=Float64[0 ,0 ,0 ,0]
    for i=1:4 
     a=real(S[i])
     b=imag(S[i])
      if b != 0
        continue
        
      end
     G[i]=a
     #println(a)
     i+= 1
    end
    Sol= filter(!iszero, G)
    =#
    indexreal=findall(x->abs.(x)<1e-12,imag.(S))

   Sol=real.(S[indexreal])
function  pickupRoot(h₁,h₂,λ₁,λ₂,Sol)
       
    for  p₁ in Sol   
         μ₁=(h₁-p₁)/(λ₁*p₁)
        
        if μ₁<=0
            continue
        else
            p₂=(λ₁*h₂)/((λ₁-λ₂)*p₁+λ₂*h₁)*p₁
            μ₂=(h₂-p₂)/(λ₂*p₂)
            if  μ₂<=0 
                continue
            end
            return   p₁,p₂
        end
    end
end 

#println(" p₁=$p₁,p₂=$p₂")

aₕ⁽ᵖ⁾=[p₁ p₂]'
aₗ⁽ᵖ⁾=Q'*( aₕ⁽ᵖ⁾- D⁻¹Qbₖ)
ηₖ⁽²⁾ =1-((aₗ⁽ᵖ⁾[1])/(normvₖ)+( uₖdotvₖ)/(vₖdotvₖ))*(aₗ⁽ᵖ⁾[2])/(normzₖ)
γₖ⁽²⁾=-1*((aₗ⁽ᵖ⁾[2])/(normzₖ))./(ηₖ⁽²⁾)
ωₖ⁽²⁾=-1*(γₖ⁽²⁾*uₖ)-vₖ
xₖ=a.+(ηₖ⁽²⁾*ωₖ⁽²⁾)

Mₖ=uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
E=1-Mₖ

#===================================================#
N=10
i = 1:1:N
i² =i.^2
a =(i².*10).+1
A=Diagonal(i²)
b=zeros(N)
α=385
x₀=sqrt(α/Quadratic(a,A,b)).*a 

#xₖ,k,E=  GenAlg(x₀,A,b,α,a,ϵ=1e-6)

#x₀=xₖ