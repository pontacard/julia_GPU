using DifferentialEquations
using LinearAlgebra
using Base.Threads



function LLG!(dS,S,params,t)
    α,B, K, γ, STTac, ω, phase ,DC_torque = params
    #println(α,B, K, γ, STTac, ω, phase ,DC_torque)
    #println(sin(S[3] + phase[3]))
    Bk = [K[1] * S[1], K[2] * S[2], K[3] * S[3]]
    #println(B+Bk+[Bac[1] * sin(ω[1] * t + phase[1]) Bac[2] * sin(ω[2] * t + phase[2]) Bac[3] * sin(ω[3] * t + phase[3])] )
    B_eff = B + Bk 
    Spin_torque = DC_torque + [STTac[1] * sin(ω[1] * t + phase[1]), STTac[2] * sin(ω[2] * t + phase[2]), STTac[3] * sin(ω[3] * t + phase[3])] 
    #println(Spin_torque)
    g = 1 / (1 + 0.288 * S[1]) #STTにつく、磁化と参照層の角度によって決まるパラメーター(磁化と参照層が平行な時に最も大きくなる、参照層が+x方向を向いたベクトルなため、磁化のx成分を見れば良い)

    
    S_cr_B = cross(S, B_eff)
    S_S_B = cross(S, S_cr_B)
    STT = g * cross(S, cross(Spin_torque, S))
    #println(STT)
    
    dS[1] = - γ * S_cr_B[1] - γ * α * S_S_B[1] - γ * STT[1]
    dS[2] = - γ * S_cr_B[2] - γ * α * S_S_B[2] - γ * STT[2]
    dS[3] = - γ * S_cr_B[3] - γ * α * S_S_B[3] - γ * STT[3]
end

function history(self, self_vec, func, tspan, S0)      #selfは辞書式、self_vecはベクトルとしてパラメーターが入ってる(DifferentialEquationsのせいでこうしてる)
    prob = ODEProblem(func, S0, tspan, self_vec)
    #println("🐰")
    result = solve(prob, dt=self.dt, adaptive=false, maxiters=1e10)
    t = result.t 
    x = [u[1] for u in result.u]
    y = [u[2] for u in result.u]
    z = [u[3] for u in result.u]

    his =transpose([x y z t])    #hisはここでしか変化させない
    #println("🐰")
    return his
end