include("../tool.jl")

using .Lya
using DifferentialEquations
using LinearAlgebra
using Base.Threads

struct paramerte
    dt
    α
    B
    BK
    γ
    ω
    phase 
end



function LLG!(dS,S,params,t)
    α,B, K, γ, SOTac, ω, phase = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))
    sinθ = sin(S[1])
    cosθ = cos(S[1])
    sinϕ = sin(S[2])
    cosϕ = cos(S[2])
    
    B_θ = cos(S[1]) * cos(S[2]) * (B[1] + K[1] * sin(S[1]) * cos(S[2])) + cos(S[1]) * sin(S[2]) * (B[2] + K[2] * sin(S[1]) * sin(S[2])) - sin(S[1]) * (B[3] + K[3] * cos(S[1]))
    #println(B_θ)
    B_ϕ = - sin(S[2]) * (B[1] + K[1] * sin(S[1]) * cos(S[2])) + cos(S[2]) * (B[2] +  K[2] * sin(S[1]) * sin(S[2]))

    SOT_θ = cosθ * cosϕ * SOTac[1] * sin(S[3] + phase[1]) + cosθ * sinϕ * SOTac[2] * sin(S[3] + phase[2]) - sinθ * SOTac[3] * sin(S[3] + phase[3])
    SOT_ϕ = - sinϕ * SOTac[1] * sin(S[3] + phase[1]) + cosϕ * SOTac[2] * sin(S[3] + phase[2])
    
    
    dS[1] = γ * B_ϕ + α * γ * B_θ + γ * SOT_θ

    dS[2] = - γ * (B_θ /sin(S[1])) + α * γ * (B_ϕ / sin(S[1])) + γ * (SOT_ϕ / sin(S[1]))
    dS[3] = ω 


end

function SOT_Bifrucation_map(paras,ax,tspan, S0,B_eval,start_step,stop_step) #step_Bは等差数列の差の値を入れる
    B_list = []
    poi_list = []
    Threads.@threads for Bac in B_eval
        #println("here") 
        his = Lya.history(paras,[0.0, Bac, 0.0],tspan,S0)       
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Poi = Lya.poincare(paras, his, ax, start_step, stop_step)
        Bac_list = fill(Bac, length(Poi))
        #print("here") 
        poi_list = append!(poi_list, Poi)
        B_list = append!(B_list, Bac_list)
        #println(Bac)
    end

    filename = "data/Bif_diagram/SOT_bifrucation_map_Bx_$(paras.B[1])_Ky_$(paras.BK[2])_$(paras.ω)GHz._start_step_$(start_step)_alpha$(α)_paper_10-35_new.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(B_list[:], poi_list[:])) # x,y,zの3列にして掃き出し
    end

end 

S0 = [1.568, 0.849, 4.833]
tspan = (0.0, 800.0)
B = [160.0, 0.0, 0.0]
BK = [0.0, 200.0, 0.0]
Bac_phase = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176335977
ω = 20.2
dt = 0.001
spin = Lya.para(dt,α,B, BK, γ, ω, Bac_phase, LLG!)

#his = Lya.history(spin,Bac,tspan,S0)
#println(length(his))
#println(Lya.poincare(spin, his, "x", 700, 790))
#Lya.Si_Sj_phase(2,1,700,750)
Bac_eval = Vector(10:0.025:35)
SOT_Bifrucation_map(spin, "y", tspan, S0, Bac_eval, 788000, 799900)