include("./../tool.jl")
using DifferentialEquations
using LinearAlgebra
using Base.Threads
using .Tool

struct STO
    dt
    α
    B
    K
    γ
    STTac
    ω
    phase 
    λ
    DC_torque
end


function LLG!(dS,S,params,t)
    #println(params)
    dt, α,B, K, γ, STTac, ω, phase, λ, DC_torque = params  #入力は交流スピントルクと直流スピントルク(SOT)
    sinθ = sin(S[1])
    cosθ = cos(S[1])
    sinϕ = sin(S[2])
    cosϕ = cos(S[2])

    g = 1 / (1 + λ * cosθ) #STTにつく、磁化と参照層の角度によって決まるパラメーター(磁化と参照層が平行な時に最も大きくなる、参照層が+x方向を向いたベクトルなため、磁化のx成分を見れば良い)
    B_θ = cosθ * cosϕ * (B[1] + K[1] * sinθ * cosϕ) +  cosθ * sinϕ * (B[2] + K[2] * sinθ * sinϕ) - sinθ * (B[3] + K[3] * cosθ)
    #println(g)
    B_ϕ = -  sinϕ * (B[1] + K[1] *  sinθ * cosϕ) + cosϕ * (B[2] + K[2] * sinθ * sinϕ)

    STT_θ = cosθ * cosϕ * (DC_torque[1] + STTac[1] * sin(S[3] + phase[1])) + cosθ * sinϕ * (DC_torque[2] + STTac[2] * sin(S[3] + phase[2])) - sinθ * (DC_torque[3] + STTac[3] * sin(S[3] + phase[3]))
    STT_ϕ = - sinϕ * (DC_torque[1] + STTac[1] * sin(S[3] + phase[1])) + cosϕ * (DC_torque[2] + STTac[2] * sin(S[3] + phase[2]))
    
    dS[1] = γ * B_ϕ + α * γ * B_θ - γ * g * STT_θ 
    dS[2] = - γ * (B_θ /sinθ) + α * γ * (B_ϕ / sinθ) - γ * g * (STT_ϕ / sinθ)
    dS[3] = ω 
    
end


function STO_Lyapunov_map(para_vec, func, per, Lya_step, cal_num, start_step, tspan, S0, STT_eval, filename) #step_Bは等差数列の差の値を入れる
    STT_list = []
    Lya_list = []
    Threads.@threads for STTac in STT_eval
        paras = STO(para_vec[1], para_vec[2], para_vec[3], para_vec[4], para_vec[5], [0, 0, STTac], para_vec[6], para_vec[7], para_vec[8], para_vec[9])
        paras_vec = [para_vec[1], para_vec[2], para_vec[3], para_vec[4], para_vec[5], [0, 0, STTac], para_vec[6], para_vec[7], para_vec[8], para_vec[9]]
        #println(paras_vec)
        Lya = Tool.matsunaga_Lyapunov(paras, paras_vec, func, per, Lya_step, cal_num, start_step, tspan, S0)
        #print(Lya) 
        Lya_list = append!(Lya_list, [Lya])
        STT_list = append!(STT_list, [STTac])
        #println(Bac)
    end

    open(filename,"w") do out
        Base.print_array(out, hcat(STT_list[:], Lya_list[:])) # x,y,zの3列にして掃き出し
    end

end 

S0 = [pi/2, 0.1, 0.0]
tspan = (0.0, 800.0)
B = [32.0, 0.0, 0.0]
BK = [0.0, 41.9, 0.0]
STTac = [0.0, 0.0, 5.0]
STTac_phase = [0.0, 0.0, 0.0]
DC_torque = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176
ω = 6.5
λ = 0.288
dt = 0.001

para_vec = [dt,α,B, BK, γ, ω, STTac_phase, λ, DC_torque]

per = [0.01,0.01, 0.01]
start_step = 700000
Lya_step = 1000
STTac_eval = Vector(1.0:0.05:10.0)

#his = Tool.history(spin, para_vec, LLG!, tspan,S0)
#a = Tool.Si_Sj_phase(his,2,1,300000,390000)
filename = "data/Lyapunov/STT_Lyapunovmap_Bx_$(B[1])_Ky_$(BK[2])_DC_torque_$(DC_torque[3])_$(ω)GHz._start_step_$(start_step)_Lyastep_$(Lya_step)_alpha$(α)_paper_1-10.txt"
STO_Lyapunov_map(para_vec, LLG!, per, Lya_step, 5, start_step, tspan, S0, STTac_eval, filename)

"""
for Bx in Bx_eval
    if BKy > Bx
        ω = γ * sqrt(BKy^2 - Bx^2)
    elseif BKy < Bx
        ω = γ * sqrt(Bx * (Bx - BKy))
    else
        Bx += 0.01
        ω = γ * sqrt(Bx * (Bx - BKy))
    end

    B = [Bx, 0.0, 0.0]
    params = paramerte(dt,α,B, BK, γ, ω, Bac_phase)
    FMR_Lyapunov_map(per,  5,params, tspan, S0, B_eval,Lya_step,start_step)
end
"""