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
       uk=A*xk+b
       vk=a-xk
       #= orthonormalize the vectors vk and uk =#
         zk=uk-((uk.*vk/ vk'vk).*vk)
         pk=vk/norm(vk)
         qk=zk/norm(zk)
         Hk=[pk,qk]
         Ak=[Hk'*A*Hk]  
         bk=Hk'*uk 
         ρk=(1/2)*[(vk'*A*vk)/(vk'*vk)+ (zk'*A*zk)/(zk'*zk)+√(((vk'*A*vk)/(vk'*vk)+(zk'*A*zk)/(zk'*zk)^2+(4*(vk'*A*zk))/(vk'*vk*zk*zk))]
         γk=(ρk)^(-1)
         ck=xk-γk*uk
         #Step3=============================================================#
         ωk=ck-a
         # calculate ηk ===========#
         ga=A*a+b
         α=q(xk)
         ηk=-((ga'*ωk)/(ω'*A*ω))-√(((ga'*ωk)/(ω'*A*ω))^2- (q(a)-α)/(ω'*A*ω))
         xkk=a+ηk*ωk
         #Step4==============================================================#
         #tolerance error condition
          Ek=1-(uk'*vk)/(norm(uk)*norm(vk))
           
          
