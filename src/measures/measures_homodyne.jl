
using LinearAlgebra, Distributions

import Base: copy
import LinearAlgebra: I,inv
using Colors, Plots

###########################

mutable struct Wigner
    n::Int64 ## number of modes
    mu::Vector{Float64}
    sigma::Matrix{Float64}
    m::Matrix{Float64}

    cst::Float64

end

function correlation_matrix(W::Wigner)
    W.m=inv(W.sigma)
    W.cst=1/(sqrt(det(2*pi*W.sigma)))
end


function Constructor(n)
    W = Wigner(n,Vector{Float64}(undef,2*n),Matrix{Float64}(undef,2*n,2*n),Matrix{Float64}(undef,2*n,2*n),1)
    
    for i in range(1,2*n)
        W.mu[i]=0
        for j in range(1,2*n)
            if(i==j)
                W.sigma[i,j]=1
            
            else
                W.sigma[i,j]=0
            end

        end
    end



    correlation_matrix(W) #compute the correlation matrix
    return W
end


function value(W::Wigner,x) #x = (x1,x2,p1,p2) compute the values

    return W.cst * exp( - transpose(x-W.mu)*W.m*(x-W.mu)  )

end

########################### CHSH









###########################

function Wigner_PseudoGauss_Integration_x(a1,a2,W1,W2) #Do the integration for sum of two gaussian wigner

    return a1*Wigner_Gaussian_Integration_x(W1)+a2*Wigner_Gaussian_Integration_x(W2)

end


function Wigner_Gaussian_Integration_x(W::Wigner) # function for (x1,x2,p1,p2)
    
    
    N = Int64 # nuber of subdivision
    a = Float64 # [-N*a,N*a] scaling

    N=40
    L=2
    a=L/N

    S = Matrix{Float64}(undef,2*N,2*N)


    for i in range(1,2*N) #initialize the matrix at 0
        for j in range(1,2*N)
            S[i,j]=0
        end

    end

    sum=0

    for i in range(-N+1,N) #integration over x
        for j in range(-N+1,N)
            for k in range(-N+1,N)
                for l in range(-N+1,N)
                    S[N+i,N+j]+=value(W,[a*i+a/2,a*j+a/2,a*k+a/2,a*l+a/2]) #center calcultion in the case of discretization
                    
                end

            end
            sum+=S[N+i,N+j] #normalize it !

        end

    end


    return S/sum

end



function integrate_p(S::Matrix{Float64},k::Int64) #integrate over pk, then return integration over the other p on [0,+infiny]

    N=40
    sum=0



    if(k!=1) 

        for i in range(1,2*N)
            for j in range(1,N-1)
                sum+=S[i,j]

            end
        end


    else

        for i in range(1,2*N)
            for j in range(1,N-1)
                sum+=S[j,i]
            end
        end




    end


return sum


end

 
W=Constructor(2)
W.mu[2]=0.5

S = Wigner_Gaussian_Integration_x(W)

using PyPlot


surf(S)
a = integrate_p(S,2)
print(a)

#pygui(true)