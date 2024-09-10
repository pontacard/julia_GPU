using DifferentialEquations
using Plots
using LinearAlgebra

struct paramerters
    dt
    σ
    ρ
    β
end


function rolenz!(dS,S,params,t)
    σ,ρ,β = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))  

    Jacobian = [-σ σ 0
                ρ-S[3] -1 -S[1]
                S[2] S[1] -β]
    #println("jaco",Jacobian)
    
    dr = transpose(reshape(S[4:end], 3,3))
    S_s = Jacobian * dr
    #println("S_s",S_s)
    #println(dr)
    #println(S_s)
    per = reshape(transpose(S_s), 1, :)
    #println("koko", per)

    dS[1] = σ * (S[2] - S[1])
    dS[2] = S[1] * (ρ - S[3]) - S[2] 
    dS[3] = S[1] * S[2] - β * S[3] 
    dS[4] = per[1]
    dS[5] = per[2]
    dS[6] = per[3]
    dS[7] = per[4]
    dS[8] = per[5]
    dS[9] = per[6]
    dS[10] = per[7]
    dS[11] = per[8]
    dS[12] = per[9]
end

function history(para_dif,tspan, S0)      #structのparasを引数にすれば良い
    lor_para  = [para_dif.σ, para_dif.ρ, para_dif.β]
    prob = ODEProblem(rolenz!, S0, tspan, lor_para)
    result = solve(prob, dt=para_dif.dt, adaptive=false)
    t = result.t 
    x = [u for u in result.u]
    #println(x[1])
    return x
end

function Lyapunov_spec(X0, step, cal_num, start_step,paras,tspan)
    dt = paras.dt
    Lya_dt = Int((tspan[2] / dt - start_step) ÷ step)
    #println(Lya_dt)
    ini_tspan = [tspan[1], dt * start_step]
    ans_ini = history(paras, ini_tspan, S0)
    #println("end",ans_ini[end][1:3])
    X0[1:3] = ans_ini[end][1:3] 
    pers  = transpose(reshape(X0[4:end], 3,3))
    #println("pers", pers)
    per_norm = [norm(pers[1 , :]) norm(pers[2 , :]) norm(pers[3 , :])]
    Lya = [0 0 0]
    for i in 1:step
        #println(i)
        tp = [0, cal_num * dt * Lya_dt]
        #println("X0", X0)
        ansp = history(paras, tp, X0)
        #println("ansp",ansp)
        per_p = transpose(reshape(ansp[Lya_dt][4:end], 3,3))
        #println("per_p", per_p)
        per_p1 = per_p[1 , :]
        per_p2 = per_p[2 , :]
        per_p3 = per_p[3 , :]
        α_i = [log(norm(per_p1)/per_norm[1]) log(norm(per_p2)/per_norm[2]) log(norm(per_p3)/per_norm[3])]
        Lya += α_i

        #println(per_p1)

        per1 = per_p1 / norm(per_p1) 
        per2_nonn = per_p2 - dot(per_p2, per1) * per1
        per2 = per2_nonn / norm(per2_nonn)
        per3_nonn = per_p3 - dot(per_p3, per1) * per1 - dot(per_p3, per2) * per2
        per3 = per3_nonn / norm(per3_nonn)
        #println("dot",dot(per1,per2), dot(per2,per3))

        per1 = per1 * per_norm[1]
        per2 = per2 * per_norm[2]
        per3 = per3 * per_norm[3]
        pers = [per1 per2 per3]
        #println("pers",pers)
        pers = reshape(transpose(pers), 1,:)
        #println("pers_af",pers)

        X_after = ansp[Lya_dt][1:3]

        X0[1:3] = X_after
        X0[4:end] = pers
        #println("hre")

        #println(X0)
        #println(Lya)
    end
    cal_time = Lya_dt * step * dt
    #println("reach ",dt)
    Lya_expo = Lya / cal_time
    println("Lya_expo ",Lya_expo)
    return Lya_expo
end 


S0 = [1.0, 1.0, 1.0, 0.01, 0.0, 0.0, 0, 0.01, 0, 0, 0, 0.01]
tspan = (0.0, 200.0)
σ = 10
ρ = 28
β = 8/3
dt = 0.001
params = paramerters(dt,σ,ρ,β)
prob = ODEProblem(rolenz!, S0, tspan, params)
start_step = 100000
Lya_step = 10000
Lyapunov_spec(S0,Lya_step,10,start_step,params,tspan)
#result = solve(prob,dt=dt, adaptive=false)
#result = solve(prob)

# 計算
#result = solve(prob, reltol=1e-12, abstol=1e-12)
"""
# 描画
t = result.t 
x = [u[1] for u in result.u]
y = [u[2] for u in result.u]
println(length(x))

plot(y[4000:length(x)], x[4000:length(x)], label="prey $x")
"""