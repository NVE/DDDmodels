#
#------------------------------------------------------------------------
#     Description:  Calculates the Shape and Scale parameter of Lambda for individual Layers. 
#
#     Author: Thomas Skaugen
#     Revised: 27.02.2019
#--------------------------------------------------------------------------


function Big2SmallLambda(GshInt,GscInt) 
    
# using LsqFit
# using Distributions
# using Statistics

arg = vcat([.1:.1:.9;],[.99])

g = Gamma(GshInt,GscInt)
avglam = quantile.(g,arg) # number of levels

trlam = zeros(Float64,10)  # Lowercase lambda for level
trlam[1] = avglam[1]

len = 20                   # length of UH  for comparison

avgUH = zeros(Float64,len)

for lag in 1:10 #Loop for 10 Layers
    
    sumwgt = sum(avglam[1:lag])
    wgt = zeros(Float64,lag)
    
    for i in 1:lag    
      wgt[i] = avglam[i]/sumwgt
    end
    avgUH = avglam[lag]*exp.(-avglam[lag]*(1:len)) # The UH to compare against
    
    dy = avgUH
    dx = [1:1:len;]
    
    #lag2
    if(lag == 2)
      @.model2(x, p) = (trlam[1]*exp(-trlam[1]*x)*wgt[1])+ (p*exp(-p*x)*wgt[2])
      p0 = [0.1]                      # parameter needs to be in brackets
      fit2 = curve_fit(model2, dx, dy, p0)  
      trlam[lag] = coef(fit2)[1]
    end

    #lag3
    if(lag == 3)
      @.model3(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2])+
      (p*exp(-p*(x))*wgt[3])  # for continuing a line in Julia, the plus (or any sign) needs to be at the end of previous line
      start =trlam[lag-1]
      p0 = [start]                    # parameter needs to be in brackets
      fit3 = curve_fit(model3, dx, dy, p0)
      trlam[lag] = coef(fit3)[1]
    end

    #lag4
    if(lag==4)
      @.model4(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2])+
         (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(p*exp(-p*(x))*wgt[4])
         start =trlam[lag-1]
         p0 = [start]                    # parameter needs to be in brackets
         fit4 = curve_fit(model4, dx, dy, p0)
         trlam[lag] = coef(fit4)[1]
    end

    #lag5
    if(lag==5)
      @.model5(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2]) +
         (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+(p*exp(-p*(x))*wgt[5])
         start = trlam[lag-1]
         p0 = [start]                    # parameter needs to be in brackets
         fit5 = curve_fit(model5, dx, dy, p0)
         trlam[lag] = coef(fit5)[1]
    end
    
    #lag6
    if(lag==6)
       @.model6(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2]) +
        (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+
        (trlam[5]*exp(-trlam[5]*(x))*wgt[5])+(p*exp(-p*(x))*wgt[6])
        start =trlam[lag-1]
        p0 = [start]                    # parameter needs to be in brackets
        fit6 = curve_fit(model6, dx, dy, p0)
        trlam[lag] = coef(fit6)[1] 
    end  

    #lag7
    if(lag==7)
      @.model7(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2])+
        (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+
        (trlam[5]*exp(-trlam[5]*(x))*wgt[5])+(trlam[6]*exp(-trlam[6]*(x))*wgt[6])+(p*exp(-p*(x))*wgt[7])
        start =trlam[lag-1]
        p0 = [start]                    # parameter needs to be in brackets
        fit7 = curve_fit(model7, dx, dy, p0)
        trlam[lag] = coef(fit7)[1] 
    end

    #lag8
    if(lag==8)
      @.model8(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2]) +
      (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+(trlam[5]*exp(-trlam[5]*(x))*wgt[5])+
      (trlam[6]*exp(-trlam[6]*(x))*wgt[6])+(trlam[7]*exp(-trlam[7]*(x))*wgt[7])+(p*exp(-p*(x))*wgt[8])
      start =trlam[lag-1]
      p0 = [start]                    # parameter needs to be in brackets
      fit8 = curve_fit(model8, dx, dy, p0)
      trlam[lag] = coef(fit8)[1] 
    end

    #lag9
    if(lag==9)
      @.model9(x, p) = (trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2])+
      (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+(trlam[5]*exp(-trlam[5]*(x))*wgt[5])+
      (trlam[6]*exp(-trlam[6]*(x))*wgt[6])+(trlam[7]*exp(-trlam[7]*(x))*wgt[7])+(trlam[8]*exp(-trlam[8]*(x))*wgt[8])+
      (p*exp(-p*(x))*wgt[9])
        start =trlam[lag-1]
        p0 = [start]                    # parameter needs to be in brackets
        fit9 = curve_fit(model9, dx, dy, p0)
        trlam[lag] = coef(fit9)[1] 
    end

    #lag10
    if(lag==10)
       @.model10(x, p) = (trlam[1]*exp(-trlam[1]*(x))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(x))*wgt[2]) +
       (trlam[3]*exp(-trlam[3]*(x))*wgt[3])+(trlam[4]*exp(-trlam[4]*(x))*wgt[4])+(trlam[5]*exp(-trlam[5]*(x))*wgt[5])+
       (trlam[6]*exp(-trlam[6]*(x))*wgt[6])+(trlam[7]*exp(-trlam[7]*(x))*wgt[7])+(trlam[8]*exp(-trlam[8]*(x))*wgt[8])+
       (trlam[9]*exp(-trlam[9]*(x))*wgt[9])+(p*exp(-p*(x))*wgt[10])
       start =trlam[lag-1]
       p0 = [start]                    # parameter needs to be in brackets
       fit10 = curve_fit(model10, dx, dy, p0)
       trlam[lag] = coef(fit10)[1] 
    end

end #loop for levels

#Estimate parameters for trlam (lambda)
modelgamma(x, p) = quantile.(Gamma(p[1], p[2]), x)
gammafit = curve_fit(modelgamma, arg, trlam, [GshInt*0.5, 2*GscInt],
                     autodiff=AutoFiniteDiff(fdjtype=Val(:central)))

#using Plots
#plot(quantile.(Gamma(coef(gammafit)[1],coef(gammafit)[2]),dy),dy)
#plot!(trlam,dy)

Gshape = coef(gammafit)[1]
Gscale = coef(gammafit)[2]

return Gshape, Gscale

end           
