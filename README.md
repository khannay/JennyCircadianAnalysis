# JennyCircadianAnalysis

This repo has the julia script file to generate the predicted DLMO times and other statistics based on wearable data.

It depends on the Arcascope repo's HCRSIMJul and WearableDataParser 


## Light Compactness Score



```{julia}
function getMidCompact(data::DataFrame, LightProxy::Function, model::CircadianModel; )

	σ=2.0*π/12.0 #pm 2 hours 
	μ=HCRSimJul.clock_circadian_convert(data.TimeTotal[1])
	sol=HCRSimJul.integrate_model(LightProxy, model, [0.70, μ ∓ σ, 0.0], 0.0, 7*24.0)
	σafter= std(Array(sol)[2,end])
	return σ/σafter-1.0 

end
```

The light compactness score measures the change in the uncertainty in the phase
variable after one week of entrainment under the lighting conditions. 

The initial distribution is taken as a gaussian distribution with stdev of 2
hours. 

In the absence of any light information L=0 then the compactness score would be
zero as the uncertainty will neither grow or decrease under a uniform propagation
forward. 

* C > 0: means that the light reduced the phase uncertainty (as the uncertainty
after goes to zero C goes to infinity.) Higher positive C would be a more regular schedule. 

* C < 0: means that the light exposure increased the uncertainty on the phase
