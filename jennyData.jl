

# Generate  average DLMO times for each subject in Jenny's data set

using HCRSimJul, WearableDataParser
using StatFiles
using DataFrames, DataFramesMeta, Dates, RollingFunctions
using CSV, Statistics
using MonteCarloMeasurements
using Plots

#demog=DataFrame(load("/Users/khannay/work/Research/Jenny Stats/demographic_sleep_PA_082218.sav") )

WearableDataParser.set_data_directory()


macro runAllJenny(myfunc,results)
	quote
		$(esc(results))=[]
		#WearableDataParser.set_data_directory() 
		jenny_acti_directory="$(WearableDataParser.data_directory)/actiwatch_no_dlmo/jenny_data/DATA" 
		myfiles=readdir(jenny_acti_directory)
		for f in myfiles
			if occursin(r".*.csv", f)
				outVal=$(esc(myfunc))(f)
				if !ismissing(outVal)
					push!($(esc(results)),outVal)
				end
			end
		end
	end
end 


function steps_mult_analysis() 

	function median_steps(filename="0") 
		#WearableDataParser.set_data_directory()
		data=WearableDataParser.read_jenny_data(filename)
		
		return median(data.Steps)

	end


	res=[]
	@runAllJenny(median_steps, res)

	pl=scatter(1:length(res), res)
	display(pl)
	

end 



function getMidCompact(data::DataFrame, LightProxy::Function, model::CircadianModel; )

	σ=2.0*π/12.0 #pm 2 hours 
	μ=HCRSimJul.clock_circadian_convert(data.TimeTotal[1])
	sol=HCRSimJul.integrate_model(LightProxy, model, [0.70, μ ∓ σ, 0.0], 0.0, 7*24.0)
	σafter= std(Array(sol)[2,end])
	return σ/σafter-1.0 

end


function genDLMOjenny(filename::String="0"; steps_multiple=2.0, makeplot=true)

	@info filename
	data=WearableDataParser.read_jenny_data(filename; makeplot=makeplot)

	if false
		median_steps = median(data.Steps[data.Steps .> 0])
		median_light = median(data.Lux[data.Lux .> 0])
		@show median_steps, median_light 
	end
	LightProxy=HCRSimJul.interp_flux(data.TimeTotal, steps_multiple .* data.Steps)
	model=SinglePopModel() 


	init=[0.70,HCRSimJul.clock_circadian_convert(data.TimeTotal[1]), 0.0]
	tend_loop=data.TimeTotal[1]+5*24.0
	sol=integrate_model(LightProxy, model,init, data.TimeTotal[1], tend_loop)
	h=0
	while h <= 10
		sol=integrate_model(LightProxy, model, sol.u[end], data.TimeTotal[1], tend_loop)
		h+=1
	end


	DLMO=integrate_observer(LightProxy, model, sol.u[end], data.TimeTotal[1], data.TimeTotal[end], model.DLMOObs)
	if length(DLMO)>0 
		Z=1/length(DLMO)*sum(exp.( im .* DLMO .* π/12.0))
		meanDLMO=12.0/π * angle(Z) 
		if meanDLMO < 0.0 
			meanDLMO += 24.0 
		end

		compactLightScore=getMidCompact(data, LightProxy, model)

		idnum=data.SubjectID[1]
		#@show  DLMO .% 24.0, meanDLMO
		R_DLMO=abs(Z)
		condition= (occursin("B", idnum)) ? "School" : "Summer"
		id=split(idnum, r"P|B")[1]
		return([idnum, id, condition, DLMO[end] % 24.0, compactLightScore, R_DLMO])
	else
		return missing
	end

end 



function makeJennyActogram(filename="0"; steps_multiple::AbstractFloat=2.0)

	data=WearableDataParser.read_jenny_data(filename)
	model=SinglePopModel()
	LightProxy=HCRSimJul.interp_flux(data.TimeTotal, data.Steps .* steps_multiple)

	print(first(data,10))
	init=[0.70,HCRSimJul.clock_circadian_convert(data.TimeTotal[1]), 0.0]
	DLMO=integrate_observer(LightProxy, model, init, data.TimeTotal[1], data.TimeTotal[end], model.DLMOObs)


	pl_actogram=HCRSimJul.plot_actogram(data; threshold=10.0) 
	HCRSimJul.plot_phase_marker!(pl_actogram, DLMO) 
	savefig("./Images/testActogramJenny_$(data.SubjectID[1]).png")
	display(pl_actogram) 

	return 1

end



function generate_file(;λ=50.0)
	res=[]
	g(f)=genDLMOjenny(f; steps_multiple=λ)
	@runAllJenny(g, res)

	# Now let's get this into a DataFrame alongside other significant covariates
	PatientNumber=[]
	Condition=[]
	DLMO=[]
	fn=[]
	R_DLMO=[]
	compactLightScore=[]


	for k in res
		push!(fn, k[1])
		push!(PatientNumber, k[2])
		push!(Condition, k[3])
		push!(DLMO, k[4])
		push!(compactLightScore, k[5])
		push!(R_DLMO, k[6])
	end

	df=DataFrame(filename=fn, PatientNumber=PatientNumber, Condition=Condition, DLMO=mean.(DLMO), Error_DLMO=std.(DLMO), CompactLS=compactLightScore, DLMOPhaseAmp=R_DLMO)
	CSV.write("DLMO_Jenny_Activity_$(Int(λ)).csv", df)
end

# Now use this in Rprogram to fit the actual LMM

generate_file(;λ=15.0)
