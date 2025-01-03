include("./../tool.jl")
using DifferentialEquations
using LinearAlgebra
using Base.Threads
using .Tool

struct FMR
    dt
    α
    B
    K
    γ
    Bac
    ω
    phase 
end

function LLG!(dS,S,params,t)
    #println(params)
    dt, α,B, K, γ, Bac, ω, phase = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))
    B_θ = cos(S[1]) * cos(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[1]) * sin(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2])) - sin(S[1]) * (B[3] + Bac[3] * sin(S[3] + phase[3]) + K[3] * cos(S[1]))
    #println(B_θ)
    B_ϕ = - sin(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2]))
    
    
    dS[1] = γ * B_ϕ + α * γ * B_θ

    dS[2] = - γ * (B_θ /sin(S[1])) + α * γ * (B_ϕ / sin(S[1]))
    dS[3] = ω 
    
end


function FMR_Lyapunov_map(para_vec, func, per, Lya_step, cal_num, start_step, tspan, S0, B_eval) #step_Bは等差数列の差の値を入れる
    B_list = []
    Lya_list = []
    Threads.@threads for Bac in B_eval
        paras = FMR(para_vec[1], para_vec[2], para_vec[3], para_vec[4], para_vec[5], [0, Bac, 0], para_vec[6], para_vec[7])
        #println(paras)
        paras_vec = [para_vec[1], para_vec[2], para_vec[3], para_vec[4], para_vec[5], [0, Bac, 0], para_vec[6], para_vec[7]]
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Lya = Tool.matsunaga_Lyapunov(paras, paras_vec, func, per, Lya_step, cal_num, start_step, tspan, S0)
        #print(Lya) 
        Lya_list = append!(Lya_list, [Lya])
        B_list = append!(B_list, [Bac])
        #println(Bac)
    end

    filename = "FMR_Lyapunovmap_Bx_$(para_vec[3][1])_Ky_$(para_vec[4][2])_$(para_vec[6])GHz._start_step_$(start_step)_Lyastep_$(Lya_step)_alpha$(α)_paper_0-25.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(B_list[:], Lya_list[:])) # x,y,zの3列にして掃き出し
    end

end 

S0 = [1.74, 0.7237, 0.0]
tspan = (0.0, 800.0)
B = [160.0, 0.0, 0.0]
BK = [0.0, 200.0, 0.0]

Bac_phase = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176335977
ω = 20.2
dt = 0.002 * 2 * pi / ω

para_vec = [dt, α, B, BK, γ, ω, Bac_phase]

per = [0.01,0.01, 0.01]
start_step = 700000
Lya_step = 1000
Bac_eval = Vector(12:3:13)

#his = Tool.history(spin, para_vec, LLG!, tspan,S0)
#a = Tool.Si_Sj_phase(his,2,1,300000,390000)
FMR_Lyapunov_map(para_vec, LLG!, per, Lya_step, 5, start_step, tspan, S0, Bac_eval)
