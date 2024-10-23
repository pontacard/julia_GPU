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


function STFMR!(dS,S,params,t)
    α,B, K, γ, Bac, SOTac, ω, phase = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))
    sinθ = sin(S[1])
    cosθ = cos(S[1])
    sinϕ = sin(S[2])
    cosϕ = cos(S[2])

    B_θ = cos(S[1]) * cos(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[1]) * sin(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2])) - sin(S[1]) * (B[3] + Bac[3] * sin(S[3] + phase[3]) + K[3] * cos(S[1]))
    #println(B_θ)
    B_ϕ = - sin(S[2]) * (B[1] + Bac[1] * sin(S[3] + phase[1]) + K[1] * sin(S[1]) * cos(S[2])) + cos(S[2]) * (B[2] + Bac[2] * sin(S[3] + phase[2]) + K[2] * sin(S[1]) * sin(S[2]))
    
    SOT_θ = cosθ * cosϕ * SOTac[1] * sin(S[3] + phase[1]) + cosθ * sinϕ * SOTac[2] * sin(S[3] + phase[2]) - sinθ * SOTac[3] * sin(S[3] + phase[3])
    SOT_ϕ = - sinϕ * SOTac[1] * sin(S[3] + phase[1]) + cosϕ * SOTac[2] * sin(S[3] + phase[2])
    
    dS[1] = γ * B_ϕ + α * γ * B_θ + γ * SOT_θ

    dS[2] = - γ * (B_θ /sin(S[1])) + α * γ * (B_ϕ / sin(S[1])) + γ * (SOT_ϕ / sin(S[1]))
    dS[3] = ω 
    
end

function history(self, Bac, SOTac, tspan, S0)      #structのparasを引数にすれば良い
    mag_para  = [self.α, self.B, self.BK, self.γ, Bac, SOTac, self.ω, self.phase]
    prob = ODEProblem(self.func, S0, tspan, mag_para)
    result = solve(prob, dt=self.dt, adaptive=false, maxiters=1e10)
    t = result.t 
    θ = [u[1] for u in result.u]
    ϕ = [u[2] for u in result.u]
    z = [u[3] for u in result.u]

    his =transpose([θ ϕ z t])    #hisはここでしか変化させない
    return his
end

function ST_FMR_Bifrucation_map(paras,ax,tspan, S0,j_eval,start_step,stop_step) #step_Bは等差数列の差の値を入れる
    j_list = []
    poi_list = []
    Threads.@threads for jac in j_eval
        Bac = 3.7 * jac         #実際はBac[mT] = jac[A/m] * 2.7 × 10^9[mT*m/A]
        SOTac = 0.153 * jac
        #println("here") 
        his = history(paras,[0.0, Bac, 0.0],[0.0, SOTac, 0.0],tspan,S0)       
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Poi = Lya.poincare(paras, his, ax, start_step, stop_step)
        jac_list = fill(jac, length(Poi))
        #println(Poi) 
        poi_list = append!(poi_list, Poi)
        j_list = append!(j_list, jac_list)
        #println(Bac)
    end
    
    filename = "data/Bif_diagram/ST_FMR_bifrucation_map_Bx_$(paras.B[1])_Ky_$(paras.BK[2])_$(paras.ω)GHz._start_step_$(start_step)_alpha$(α)_paper_0-2.7.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(j_list[:], poi_list[:])) # x,y,zの3列にして掃き出し
    end
    
end 

S0 = [pi/2, 0, 0.0]
tspan = (0.0, 800.0)
B = [180.0, 0.0, 0.0]
BK = [0.0, 200.0, -1000.0]
jac_phase = [0.0, 0.0, 0.0]
α = 0.01
γ = 0.176335977
ω = 37.16
dt = 0.001
spin = Lya.para(dt,α,B, BK, γ, ω, jac_phase, STFMR!)
println(ω)

#his = Lya.history(spin,Bac,tspan,S0)
#println(length(his))
#println(Lya.poincare(spin, his, "x", 700, 790))
#Lya.Si_Sj_phase(2,1,700,750)
jac_eval = Vector(0.1:0.005:2.7)
#jac_eval = [2.7]
ST_FMR_Bifrucation_map(spin, "y", tspan, S0, jac_eval, 788000, 799900)
