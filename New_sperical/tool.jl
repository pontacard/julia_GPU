
module Tool
using DifferentialEquations
using LinearAlgebra
using Plots
using DelimitedFiles

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

function Si_Sj_phase(his,i,j,start_step, end_step)
    x = his[i,:]
    y = his[j,:]
    println(his[:,1])
    display(plot(x[start_step:end_step], y[start_step:end_step]))
end

function history_text(his,filename)
    open(filename, "w") do file
        writedlm(file,his)
    end
end

function poincare(self, his, ax, start_step, end_step)
    S = []
    Sx(θ,φ) = sin(θ) * cos(φ) 
    Sy(θ,φ) = sin(θ) * sin(φ) 
    Sz(θ) = cos(θ)
    if ax == "x"
        S = Sx.(his[1,:], his[2,:])
    elseif ax == "y"
        S = Sy.(his[1,:], his[2,:])
    elseif ax == "z"
        S = Sz.(his[1,:])
    else
        error(DomainError, ": ax is considering only x or y or z")
     #ユークリッド空間でのスピンの軌跡を格納完了
    end

    dt = self.dt
    T = 2 * pi / self.ω
    storobo = Int(T ÷ dt)
    println(storobo)
    println(typeof(start_step))

    return S[start_step:storobo:end_step]
end

function x_y_poincare(self, his, start_step, end_step, filename = "")
    dt = self.dt
    T = 2 * pi / self.ω
    storobo = Int(T ÷ dt) + 1
    println(storobo)
    x_list = his[1,:]
    y_list = his[2,:]
    poi_x = x_list[start_step:storobo:end_step]
    poi_y = y_list[start_step:storobo:end_step]
    list = transpose([poi_x poi_y])
    if filename != ""
        open(filename, "w") do file
            println("🍓")
            writedlm(file,list)
        end
    end
    display(plot(poi_x, poi_y, st=:scatter, mc=:blue, ms=1.3, ma=0.5, msw = 0))
end

function FMR_Lyapunov_map(self, per, cal_num,tspan, S0,sta_B,end_B,step_B,Lya_step,start_step) #step_Bは等差数列の差の値を入れる
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
    

function the_phi_poincare(self, his, start_step, end_step, savefig = false, filename = "")
    dt = self.dt
    T = 2 * pi / self.ω
    storobo = Int(T ÷ dt)
    println(storobo)
    the_list = his[1,:]
    phi_list = his[2,:]
    poi_the = the_list[start_step:storobo:end_step]
    poi_phi = phi_list[start_step:storobo:end_step]
    list = transpose([poi_the poi_phi])
    if savefig
        open(filename, "w") do file
            writedlm(file,list)
        end
    end
    display(plot(poi_phi, poi_the, st=:scatter, mc=:blue, ms=2, ma=0.5, msw = 0))
end

function cos_phi_poincare(self, his, start_step, end_step)
    Sz(θ) = cos(θ)
    S = Sz.(his[1,:])

    dt = self.dt
    T = 2 * pi / self.ω
    storobo = Int(T ÷ dt)
    println(storobo)
    phi_list = his[2,:]
    poi_cos = S[start_step:storobo:end_step]
    poi_phi = phi_list[start_step:storobo:end_step]
    display(plot(poi_phi, poi_cos, st=:scatter, mc=:blue, ms=1.3, ma=0.5, msw = 0))
end

function FMR_Lyapunov_map(self, per, cal_num,tspan, S0,sta_B,end_B,step_B,Lya_step,start_step) #step_Bは等差数列の差の値を入れる
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
   
function matsunaga_Lyapunov(self, self_vec, func, pertu, step, cal_num, start_step, tspan, S0)
    dt = self.dt
    Lya_dt = Int((tspan[2] / dt - start_step) ÷ step)
    #println(Lya_dt)
    ans0 = history(self, self_vec, func, tspan, S0)
    #println(ans0[: , start_step])
    

    dX0 = ans0[: , start_step][1:3] + pertu
    #println(dX0)
    ansp = history(self, self_vec, func, [0, cal_num * dt * Lya_dt], dX0)
    #println(ansp[:, 1])

    dist0 = norm(ans0[: , start_step][1:3] - ansp[:, 1][1:3])
    #println(ans0[: , start_step] - ansp[:, 1], dist0)

    Lya = 0.0
    for i in 1:step
        end_st = Int(i * Lya_dt + start_step)
        #ansp[:, Lya_dt + 1][3] = ans0[: , end_st][3]
        #println(ansp[:, Lya_dt + 1]," ", ans0[: , end_st])
        p_i = norm(ans0[: , end_st][1:3] - ansp[:, Lya_dt + 1][1:3]) / dist0
        #println("here")
        per_X0i = ans0[: , end_st][1:3] + (ansp[:, Lya_dt + 1][1:3] - ans0[: , end_st][1:3]) / p_i
        tp = [0, cal_num * dt * Lya_dt]
        ansp = history(self, self_vec, func, tp, per_X0i)

        Lya += log(p_i)
        #println(Lya)
    end

    cal_time = Lya_dt * step * dt
    Lya_expo = Lya / cal_time
    return Lya_expo
end 

end