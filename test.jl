using Base.Threads

"""
iow = open("julia_GPU/data/te.txt","w")
write(iow,string([1 2]))
close(iow)

open("julia_GPU/data/SOT_Lyapunovmap_Bx_100.0_Ky_200_30.542287136629696GHz._start_step_70000000_Lyastep_1001_paper10-25.txt","r") do ior
    text = read(ior,String)
    println(text)
    text = replace(text,"\n" => "") #改行コードの除去
    global val = text
    val2 = text
    #println(val2)
  end

x = [1 2 4 5 6]
y = [3 4 2 1 4]
open("julia_GPU/data/file.txt","w") do out
    Base.print_array(out, hcat(x[:],y[:])) # x,y,zの3列にして掃き出し
end
"""

a = 1
s = "ass$a"
println(s)

@show nthreads()



  