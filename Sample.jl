using LinearAlgebra
using Statistics
#Algorithm 1:=  Maximal 2-dimensional inside ball algorithm.
#= 1. Given x^0\in Ω(ξ), ϵ >0, and set k=0.
   2.  Calculate  \gamma some how ( * below)! and \u_k  and c_k.
   3.Calculate ω_k and η_k, then x^{k+1}.
   4. if tolerence condition (2.8) does NOT holds, set k=k+1 and go to step 2. =#
   #=============================================================#
   #= formulas:
    -ξ :={x\in \R^{10} | Σ_{1}^{10} i^2x_{i}^{2}=385}
    - u_k=∇q(x_k)/2 = Ax_k + b.
    - c_k = x_k − γ_k*u_k,

   =#
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
      x0=sqrt(385/q(a)).*a #q(x0)=384.99999999999994, to be sure that inital point belongs to Ω#
      e=(10)^(-6) 

      #Step2 ==========================================================================#
       u0=A*x0
       v0=a-x0

       j=1:1000
       for k in j
            uk=A*xk
            vk=a-xk
            Mk=(uk'*vk)/(norm(uk)* norm(vk))  #Mk is qoutient in tolerance condition
            E1=1-Mk[1]                        #tolerance condition
           if E1<=e
              
               println(" The algorithm find the projection onto ξ after n=$k itteration, x=$xk  ")
            break
           end 
      #= orthonormalize the vectors vk and uk =#
             zk=uk-((uk'*vk)/ (vk'*vk)*vk)
             pk=vk/norm(vk)
             qk=zk/norm(zk)
             Hk=[pk qk]   #H \in R^{n*2}
             Ak=[Hk'*A*Hk]  #Ak \in R^{2*2b}
             bk=Hk'*uk      # bk \in 2*1 vector
     #  spectral radius========================
              M1=(vk'*A*vk)/(vk'*vk)
              M2=(zk'*A*zk)/(zk'*zk)
              M3=((M1-M2)^2 )
              M4=4*(vk'*A*zk)^2/(vk'*vk*zk'*zk)
              M5=sqrt(M3+M4)
              pk=[M1+M2+M3]/2    #spectral radius (ρ)
              tk=(1)/(pk)    # inverse of spectral radius (γ)
              ck=xk-(tk.*uk)
     # xk+1   Step3  =============================
              wk=ck-a
              # ηk =============================
              ga=A*a
              o=q(x0)    #α = q(xk+1 ) ????
              N1=(ga'*wk)/(wk'*A*wk)
              N2=(N1).^2
              N3=(q(a)-o)./(wk'*A*wk)
              hk=-N1-sqrt(N2-N3)  #η is a number too close to 1.
              xk=a+(hk.*wk)
           println("Itt :$k,  x=$xk" )
       end
                
         
         
          
         
    