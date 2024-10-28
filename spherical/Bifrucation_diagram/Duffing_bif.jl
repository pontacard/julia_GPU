using DifferentialEquations
using LinearAlgebra
using Base.Threads
using Plots

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
    dX[2] = - α * X[2] + β * X[1] - γ * X[1] * X[1] * X[1] + fac * cos(X[3]) 
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

    his =transpose([x y z t])    #hisはここでしか変化させない
    return his
end

function poincare(self, his, ax, start_step, end_step)
    if ax == "x"
        S = his[1,:]
    elseif ax == "y"
        S = his[2,:]
    else
        error(DomainError, ": ax is considering only x or y")
     #ユークリッド空間でのスピンの軌跡を格納完了
    end

    dt = self.dt
    T = 2 * pi / self.ω
    storobo = Int(T ÷ dt)
    #println(storobo)
    #println(typeof(start_step))

    return S[start_step:storobo:end_step]
end

function Duf_Bifrucation_map(paras,ax,tspan, S0,f_eval,start_step,stop_step) #step_Bは等差数列の差の値を入れる
    f_list = []
    poi_list = []
    Threads.@threads for fac in f_eval
        #println("here") 
        his = history(paras,fac,tspan,S0)       
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Poi = poincare(paras, his, ax, start_step, stop_step)
        println(Poi)
        fac_list = fill(fac, length(Poi))
        #print("here") 
        poi_list = append!(poi_list, Poi)
        f_list = append!(f_list, fac_list)
        #println(Bac)
    end

    filename = "data/Bif_diagram/Duf_bifrucation_map_α_$(paras.α)_β_$(paras.β)_γ_$(paras.γ)_$(paras.ω)GHz._start_step_$(start_step)_paper_0-25.txt"
    open(filename,"w") do out
        Base.print_array(out, hcat(f_list[:], poi_list[:])) # x,y,zの3列にして掃き出し
    end

end 

X0 = [0.859, 0.0, 0.0]
tspan = (0.0, 800.0)
α = 1
β = 130
γ = 176
ω = 16.11
dt = 0.001
fac = 30.1
phase = 0.0

spin = paramerte(dt,α, β, γ, ω, phase, Duffing!)
his = history(spin, fac, tspan, X0)

#println(his[1])
#x = his[1,:]
#y = his[2,:]
#plot(x[700000:750000], y[700000:750000])

#println(ω)
fac_eval = Vector(30:0.1:80)
Duf_Bifrucation_map(spin, "x", tspan, X0, fac_eval, 788000, 799900)
#plot(his[1][900000:950000], his[2][900000:950000])

