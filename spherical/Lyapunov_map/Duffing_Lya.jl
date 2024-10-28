using DifferentialEquations
using LinearAlgebra
using Base.Threads

struct paramerte
    dt
    α
    β
    γ
    ω
    phase 
    func
end


function Duffing!(dX,X,params,t)
    α, β, γ, fac, ω, phase = params
    #println(sin(S[3] + phase[1]))
    #println(sin(S[3] + phase[3]))
    
    dX[1] = X[2]
    dX[2] = - α * X[2] + β * X[1] - γ * X[1]^3 + fac * sin(X[3]) 
    dX[3] = ω 
    
end

function history(self, fac, tspan, S0)      #structのparasを引数にすれば良い
    mag_para  = [self.α, self.β, self.γ, fac, self.ω, self.phase]
    prob = ODEProblem(self.func, S0, tspan, mag_para)
    result = solve(prob, dt=self.dt, adaptive=false, maxiters=1e10)
    t = result.t 
    x = [u[1] for u in result.u]
    y = [u[2] for u in result.u]
    z = [u[3] for u in result.u]

    his =[x y z]  #hisはここでしか変化させない
    return his
end


function Duf_Lyapunov_map(per, cal_num,paras,tspan, X0,f_eval,Lya_step,start_step) #step_Bは等差数列の差の値を入れる
    f_list = []
    Lya_list = []
    Threads.@threads for fac in f_eval
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Lya = matsunaga_Lyapunov(per, Lya_step, cal_num, start_step, paras, fac, tspan, X0)
        print(Lya) 
        Lya_list = append!(Lya_list, [Lya])
        f_list = append!(f_list, [fac])
        #println(fac)
    end

    filename = "data/Lyapunov/Duf_Lyapunovmap_α_$(paras.α)_β_$(paras.β)_γ_$(paras.γ)_$(paras.ω)GHz._start_step_$(start_step)_Lyastep_$(Lya_step)_paper.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(f_list[:], Lya_list[:])) # x,y,zの3列にして掃き出し
    end

end 
    
function matsunaga_Lyapunov(pertu, step, cal_num, start_step,paras,fac, tspan, X0)
    dt = paras.dt
    Lya_dt = Int((tspan[2] / dt - start_step) ÷ step)
    #println(Lya_dt)
    ans0_ = history(paras,fac, tspan, X0)
    ans0 = transpose(ans0_)
    println(ans0[: , start_step])
    

    dX0 = ans0[: , start_step] + pertu
    #println(dX0)
    ansp = transpose(history(paras, fac, [0, cal_num * dt * Lya_dt], dX0))
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
        ansp = transpose(history(paras, fac, tp, per_X0i))

        Lya += log(p_i)
        #println(Lya)
    end

    cal_time = Lya_dt * step * dt
    Lya_expo = Lya / cal_time
    #println(Lya_expo)
    return Lya_expo
end 

X0 = [0.6, 0.0, 0.0]
tspan = (0.0, 800.0)
α = 1
β = 130
γ = 176
ω = 16.11
dt = 0.001
fac = 38.6
phase = 0.0
per = [0.03, 0.03, 0.03]
start_step = 700000
Lya_step = 1001
fac_eval = Vector(30:0.1:80)

params = paramerte(dt,α,β, γ, ω, phase, Duffing!)
Duf_Lyapunov_map(per,  5,params, tspan, X0, fac_eval,Lya_step,start_step)