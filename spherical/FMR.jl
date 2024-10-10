include("./tool.jl")

using .Lya
using DifferentialEquations
using LinearAlgebra
using Plots

function LLG!(dS,S,params,t)
    α,B, K, γ, Bac, ω, phase = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))
    B_θ = cos(S[1]) * cos(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[1]) * sin(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2])) - sin(S[1]) * (B[3] + Bac[3] * sin(S[3] + phase[3]) + K[3] * cos(S[1]))
    #println(B_θ)
    B_ϕ = - sin(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2]))
    
    
    dS[1] = γ * B_ϕ + α * γ * B_θ

    dS[2] = - γ * (B_θ /sin(S[1])) + α * γ * (B_ϕ / sin(S[1]))
    dS[3] = ω 
    
end

S0 = [pi/2, 0.0, 0.0]
tspan = (0.0, 800.0)
B = [160.0, 0.0, 0.0]
BK = [0.0, 200.0, 0.0]
Bac = [0.0, 10.92, 0.0]
Bac_phase = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176335977
ω = 21.16
dt = 0.001
spin = Lya.para(dt,α,B, BK, γ, ω, Bac_phase, LLG!)
per = [0.01,0.01, 0.01]
start_step = 700000
Lya_step = 1000

path = "/Users/tatsumiryou/PycharmProjects/spin/others/plot_from_julia/data/phase/"
filename = "test.txt"
Lya.history(spin,Bac,tspan,S0)
Lya.Si_Sj_phase(2,1,700000,750000)
Lya.history_text(path * filename)

#@time Lya.matsunaga_Lyapunov(spin, per, Lya_step, 5, start_step, Bac, tspan, S0)