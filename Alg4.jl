m1=1
m2=1
c1=0.1
c2=0.8
max_itr=20
k=1
#i=0 
while k<=max_itr
     b=mod(k,c1+c2) # Remainder after division (modulo operation)
        if b<c1
            println("use Alg 2 to find γₖ⁽²⁾")
        else 
            println(" γₖ=c1*γₖ⁽¹⁾+ c2*γₖ⁽²⁾")
        end
 k+=1
continue  
end