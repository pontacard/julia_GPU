
module Lya
using DifferentialEquations
using LinearAlgebra
using Plots

 struct para
    dt
    α
    B
    BK
    γ
    ω
    phase 
    func
 end

function history(self, Bac, tspan, S0)      #structのparasを引数にすれば良い
    mag_para  = [self.α, self.B, self.BK, self.γ, Bac, self.ω, self.phase]
    prob = ODEProblem(self.func, S0, tspan, mag_para)
    result = solve(prob, dt=self.dt, adaptive=false)
    t = result.t 
    θ = [u[1] for u in result.u]
    ϕ = [u[2] for u in result.u]
    z = [u[3] for u in result.u]
    return [θ ϕ z]
end

function FMR_Lyapunov_map(self,per, cal_num,tspan, S0,sta_B,end_B,step_B,Lya_step,start_step) #step_Bは等差数列の差の値を入れる
    B_eval = Vector(sta_B :step_B: end_B)
    B_list = []
    Lya_list = []
    for Bac in B_eval
        #duf = FMR(t_span,α, B_ex,BK, γ,[0,B,0], ω,phase, S0)
        Lya = matsunaga_Lyapunov(self, per, Lya_step, cal_num, start_step, [0.0, Bac, 0.0], tspan, S0)
        #print(Lya) 
        Lya_list = append(Lya_list, [Lya])
        B_list = append(B_list, [Bac])
    end
    #np.savetxt(f"csv/maps_0.01_light/FMR_Lyapunovmap_Bx_{B_ex[0]}_Ky_{K[1]}_{omega}GHz._start_step_{start_step}_Lyastep_{Lya_step}_alpha{alpha}_paper_0-25.txt", np.stack([B_list,Lya_list]))
end 
    
function matsunaga_Lyapunov(self,pertu, step, cal_num, start_step,Bac, tspan, S0)
    dt = self.dt
    Lya_dt = Int((tspan[2] / dt - start_step) ÷ step)
    println(Lya_dt)
    ans0_ = history(self,Bac, tspan, S0)
    ans0 = transpose(ans0_)
    println(ans0[: , start_step])
    

    dX0 = ans0[: , start_step] + pertu
    #println(dX0)
    ansp = transpose(history(self, Bac, [0, cal_num * dt * Lya_dt], dX0))
    println(ansp[:, 1])

    dist0 = norm(ans0[: , start_step][1:2] - ansp[:, 1][1:2])
    println(ans0[: , start_step] - ansp[:, 1], dist0)

    Lya = 0.0
    for i in 1:step
        end_st = Int(i * Lya_dt + start_step)
        #ansp[:, Lya_dt + 1][3] = ans0[: , end_st][3]
        println(ansp[:, Lya_dt + 1]," ", ans0[: , end_st])
        p_i = norm(ans0[: , end_st] - ansp[:, Lya_dt + 1]) / dist0
        #println("here")
        per_X0i = ans0[: , end_st] + (ansp[:, Lya_dt + 1] - ans0[: , end_st]) / p_i
        tp = [0, cal_num * dt * Lya_dt]
        ansp = transpose(history(self, Bac, tp, per_X0i))

        Lya += log(p_i)
        #println(Lya)
    end

    cal_time = Lya_dt * step * dt
    Lya_expo = Lya / cal_time
    return Lya_expo
end 

end
