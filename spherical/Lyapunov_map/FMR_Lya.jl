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

function history(para_dif, Bac, tspan, S0)      #structのparasを引数にすれば良い
    mag_para  = [para_dif.α, para_dif.B, para_dif.BK, para_dif.γ, Bac, para_dif.ω, para_dif.phase]
    prob = ODEProblem(LLG!, S0, tspan, mag_para)
    result = solve(prob, dt=para_dif.dt, adaptive=false)
    t = result.t 
    θ = [u[1] for u in result.u]
    ϕ = [u[2] for u in result.u]
    z = [u[3] for u in result.u]
    return [θ ϕ z]
end

function FMR_Lyapunov_map(per, cal_num,paras,tspan, S0,B_eval,Lya_step,start_step) #step_Bは等差数列の差の値を入れる
    B_list = []
    Lya_list = []
    Threads.@threads for Bac in B_eval
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Lya = matsunaga_Lyapunov(per, Lya_step, cal_num, start_step, paras, [0.0, Bac, 0.0], tspan, S0)
        #print(Lya) 
        Lya_list = append!(Lya_list, [Lya])
        B_list = append!(B_list, [Bac])
        #println(Bac)
    end

    filename = "data/Lyapunov/FMR_Lyapunovmap_Bx_$(paras.B[1])_Ky_$(paras.BK[2])_$(paras.ω)GHz._start_step_$(start_step)_Lyastep_$(Lya_step)_alpha$(α)_paper_0-25.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(B_list[:], Lya_list[:])) # x,y,zの3列にして掃き出し
    end

end 
    
function matsunaga_Lyapunov(pertu, step, cal_num, start_step,paras,Bac, tspan, S0)
    dt = paras.dt
    Lya_dt = Int((tspan[2] / dt - start_step) ÷ step)
    #println(Lya_dt)
    ans0_ = history(paras,Bac, tspan, S0)
    ans0 = transpose(ans0_)
    #println(ans0[: , start_step])
    

    dX0 = ans0[: , start_step] + pertu
    #println(dX0)
    ansp = transpose(history(paras, Bac, [0, cal_num * dt * Lya_dt], dX0))
    #println(ansp[:, 1])

    dist0 = norm(ans0[: , start_step] - ansp[:, 1])
    #println(ans0[: , start_step] - ansp[:, 1], dist0)

    Lya = 0.0
    for i in 1:step
        end_st = Int(i * Lya_dt + start_step)
        #ansp[:, Lya_dt + 1][3] = ans0[: , end_st][3]
        #println(ansp[:, Lya_dt + 1]," ", ans0[: , end_st])
        p_i = norm(ans0[: , end_st] - ansp[:, Lya_dt + 1]) / dist0
        #println("here")
        per_X0i = ans0[: , end_st] + (ansp[:, Lya_dt + 1] - ans0[: , end_st]) / p_i
        tp = [0, cal_num * dt * Lya_dt]
        ansp = transpose(history(paras, Bac, tp, per_X0i))

        Lya += log(p_i)
        #println(Lya)
    end

    cal_time = Lya_dt * step * dt
    Lya_expo = Lya / cal_time
    #println(Lya_expo)
    return Lya_expo
end 


# 問題定義 (地球・月系のNear-Rectilinear Halo Orbitを与える初期条件)
S0 = [pi/2, 0.0, 0.0]
tspan = (0.0, 800.0)
BK = [0.0, 0.0, 200.0]
Bac = [0.0, 15, 0.0]
Bac_phase = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176335977
dt = 0.001
per = [0.01, 0.01, 0.01]
start_step = 700000
Lya_step = 1001
B_eval = Vector(0:0.2:60)

Bx_eval = Vector(10:1:205)
BKy = BK[2]

ω = 35.2
B = [0.0, 0.0, 0.0]
params = paramerte(dt,α,B, BK, γ, ω, Bac_phase)
FMR_Lyapunov_map(per,  5,params, tspan, S0, B_eval,Lya_step,start_step)

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