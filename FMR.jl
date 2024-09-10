struct FMR
    tspan
    α
    B
    BK
    γ
    Bac
    ω
    phase 
    S0 
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


using DifferentialEquations


# 問題定義 (地球・月系のNear-Rectilinear Halo Orbitを与える初期条件)
S0 = [pi/2, 0.0, 0.0]
tspan = (0.0, 1000000.0)
B = [160.0, 0.0, 0.0]
BK = [0.0, 200.0, 0.0]
Bac = [0.0, 5, 0.0]
Bac_phase = [0.0, 0.0, 0.0]
α = 0.05
γ = 0.176335977
ω = 20.232
dt = 0.03
params = [α,B, BK, γ, Bac, ω, Bac_phase]
prob = ODEProblem(LLG!, S0, tspan, params)
@time result = solve(prob,dt=dt, adaptive=false)
#result = solve(prob)

# 計算
#result = solve(prob, reltol=1e-12, abstol=1e-12)

println(length(x))



#plot(y[140000:length(x)], x[140000:length(x)], label="prey $x")
