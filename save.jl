using LinearAlgebra
using PolynomialRoots
function Quadratic(x,A,b)
    return dot(x,A*x)+2*dot(b,x)
end


 function GenAlg(x₀,A,b,α,a ;ϵ=1e-6, max_itr=1000 )
    xₖ=x₀
    gₐ=A*a
    qa_α = Quadratic(a,A,b)-α
    uₖ =A*xₖ+b
    vₖ =a-xₖ 
    uₖdotvₖ=dot(uₖ,vₖ)
    vₖdotvₖ=dot(vₖ,vₖ)
    normuₖ=norm(uₖ)
    normvₖ=sqrt(vₖdotvₖ)
    Mₖ =uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
    E=1-Mₖ 
    k=0
    while k<= max_itr && E>=ϵ
       
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
        M1=(vₖdotAvₖ)/(vₖdotvₖ)
        M2=(zₖdotAzₖ)/(zₖdotzₖ)
        M3=((M1-M2)^2 )
        M4=4*(vₖdotAzₖ)^2/(vₖdotvₖ*zₖdotzₖ)
        M5=sqrt(M3+M4)
        ρₖ=[M1+M2+M5]/2 
        γₖ=1/ρₖ
        cₖ=xₖ-( γₖ.*uₖ)  #the center of maximal 2-d inside ball
        #  Step3  ==== Calculation of  xk+1  ==============================
        wₖ=cₖ-a
        # ηk =============================
        Awₖ=A*wₖ
        wₖdotAwₖ= dot(wₖ,Awₖ)

         
        N1=dot(gₐ,wₖ)/(wₖdotAwₖ)
        N2=(N1).^2
        N3=(qa_α)./(wₖdotAwₖ)
        ηₖ=-N1-sqrt(N2-N3)  #η is a number too close to 1.
        xₖ=a+(ηₖ.*wₖ)  
        uₖ =A*xₖ+b
        vₖ =a-xₖ 
        uₖdotvₖ=dot(uₖ,vₖ)
        vₖdotvₖ=dot(vₖ,vₖ)
        normuₖ=norm(uₖ)
        normvₖ=sqrt(vₖdotvₖ)
        Mₖ =uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
        E=1-Mₖ 
        k+=1
       #Alg2================================================
        aₗ=[normvₖ,0]
        Aₖ=((vₖdotAvₖ)/ (normvₖ^2),(vₖdotAzₖ)/(normvₖ*normzₖ);(zₖdotAvₖ)/(normvₖ*normzₖ),(zₖdotAzₖ)/(normzₖ^2)   )
        bₖ=((uₖdotvₖ)/(normvₖ);(uₖdotzₖ)/(normzₖ))
        D, Q = eigen(A)
        Qbₖ=Q*bₖ
        D⁻¹Qbₖ=D.\Qbₖ
        #Aₖ=Q'*D*Q
        β=dot(Qbₖ,D⁻¹Qbₖ)
        aₕ=Q*aₗ+D⁻¹Qbₖ
        #Solve Quadratic=========================================
        h₁=aₕ[1]
        h₂=aₕ[2]
        λ₁=D[1,1]
        λ₂=D[2,2]
        β=2
        #solve Quadratic in lemma 1 to find p₁ then p₂ with following formula
        P=((λ₁*(λ₁-λ₂)^2))*(p₁)^4+(2*λ₁*λ₂*h₁*(λ₁-λ₂))*(p₁)^3+((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))*(p₁^2)-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)*(p₁)-(β*λ₂^2*h₁^2)
            a=((λ₁*(λ₁-λ₂)^2))
            b=(2*λ₁*λ₂*h₁*(λ₁-λ₂))
            c=((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))
            d=-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)
            f=-(β*λ₂^2*h₁^2)
        S=roots([a, b, c, d, f])
        
        #choose the right root
        function  pickupRoot(h₁,λ₁,S)
            for p₁ in S                 
                 μ₁=(h₁-p₁)/(λ₁*p₁)
                
                if μ₁<=0
                    continue
                else
                    p₂=(λ₁*h₂)/((λ₁-λ₂)*p₁+λ₂*h₁)
                    μ₂=(h₂-p₂)/(λ₂*p₂)
                    if μ₁ !≈ μ₂

                        continue
                    end
                    return   p₁,p₂
                end
            end
           
        end

        
        aₗ⁼=Q'( aₕ⁼- D⁻¹Qbₖ)
       # xₖ₊₁ ===========================================================#
        ηₖ=1-((a⁼ₗ[1])/(normvₖ)+( uₖdotvₖ)/(vₖdotvₖ))*(a⁼ₗ[2])/(normzₖ)
        γₖ=-(a⁼ₗ[2])/(normzₖ)*(1/ηₖ)

    end
    return xₖ,k,E
end
     
     N=10
     i = 1:1:N
     i² =i.^2
     a =(i².*10).+1
     A=Diagonal(i²)
     b=zeros(N)
     α=385
     x₀=sqrt(α/Quadratic(a,A,b)).*a 
    
     xₖ,k,E=  GenAlg(x₀,A,b,α,a,ϵ=1e-6)



     #=Alg3
     γ₃ₖ₊ᵢ=γ₃ₖ₊ᵢ⁽²⁾     for i=0,1
     γ₃ₖ₊₂=γ₃ₖ₊ᵢ⁽¹⁾ 
      =#

      #=Alg4
      m₁ =1
      m₂ =1
      c₁ =0.1
      c₂ =0.8

      while k<= max_itr && mod(k,m₁+m₂<=m₁)
        γₖ=γₖ⁽²⁾
      else
        γₖ=c₁*γₖ⁽¹⁾+c₂*γₖ⁽²⁾
      end

     =#