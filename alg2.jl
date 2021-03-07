using LinearAlgebra,plots,random, Statistics

#Algorithm 2:=  Sequential 2-dimensional projection algorithm.
#= 1. Given x^0\in Ω(ξ), ϵ >0, and set k=0.
   2. calculate a_l^* as projection of a onto two dimensional ellipsoid ξ_k
   3.Calculate \gamma_k and η_k, then x^{k+1}.
   4. if tolerence condition (2.8) does NOT holds, set k=k+1 and go to step 2. =#
   #=============================================================#
   
   #= Ask Inputs=#
   println(" Please inter the Tolerance criterion (ϵ):")

   #= variables============================================================#
     i = 1:1:10
     i2 =i.^2
     a =(i2.*10).+1
     A=Diagonal(i2)
     b=0
#Ellipsoid formula======================================================#
     x = 1:1:10 #Limit Point#
     function q(x)
       ret = sum((i.^2).*(x.^2))
      return ret 
       end

      
      #Step1 ==========================================================================#
      k=0 
      x₀ =sqrt(385/q(a)).*a #q(x0)=384.99999999999994, to be sure that inital point belongs to Ω#
      e=(10)^(-6) 

      #Step2 ==========================================================================#
       u₀ =A*x₀
       v₀ =a-x₀

       j=1:1000
       for k in j
            uₖ =A*xₖ
            vₖ =a-xₖ
            Mₖ =(uₖ'*vₖ)/(norm(uₖ)* norm(vₖ))  #Mk is qoutient in tolerance condition
            E=1-Mₖ[1]                        #tolerance condition
           if E<=e
              
               println(" The algorithm find the projection onto ξ after n=$k iteration, x=$xk  ")
            break
           end 
      #= orthonormalize the vectors vk and uk =#
             zₖ =uₖ-((uₖ'*vₖ)/ (vₖ'*vₖ)*vₖ)
             pₖ =vₖ/norm(vₖ)
             qₖ =zₖ/norm(zₖ)
             Hₖ =[pₖ qₖ]   #H \in R^{n*2}
             Aₖ=[Hₖ'*A*Hₖ]  #Ak \in R^{2*2b}
             bₖ=Hₖ'*uₖ     # bk \in 2*1 vector
             #compute projection of aₗ onto 2-dimensional ellipsoid ξₖˡ====================
               
            #  QDQ eigendecomposition,
             val, vec = eigen(Aₖ)    #val, vec = eigen(A)
             D= Diagonal(val)
             Q=vec
             Aₖ=Q'*D*Q

              aₗ=(norm(vₖ),0)'# a corresponds to aₗ in in l-space
              aₕ=Q*aₗ+(inv(D)*Q*bₖ)
              aₕᵖ=(h₁ᵖ,h₂ᵖ)  # solve the quadratic in lemma 1 to get a_h^* then ==>
              aₗᵖ=Q'( aₕᵖ-inv(D)*Q*bₖ)

     # xk+1   Step3  =============================
              
              # ηk =============================
             ηₖ=1-(h₁ᵖ)/(norm(vₖ))+(uₖ'*vₖ)/(vₖ'*vₖ)*(h₂ᵖ)/(norm(zₖ))
             γₖ=-1*(h₂ᵖ)/(norm(zₖ))*(1)/(ηₖ)
             xₖ=xₖ-ηₖ*γₖ*uₖ+(1-ηₖ)vₖ  #x\_k+1
           println("It :$k,  x=$xk" )
       end
                
         
         
          
         
    