

# Generate  average DLMO times for each subject in Jenny's data set

using HCRSimJul, WearableDataParser
using StatFiles
using DataFrames, DataFramesMeta, Dates, RollingFunctions
using CSV, Statistics
using MonteCarloMeasurements
using Plots

#demog=DataFrame(load("/Users/khannay/work/Research/Jenny Stats/demographic_sleep_PA_082218.sav") )




macro runAllJenny(myfunc,results)
	quote
		$(esc(results))=[]
		WearableDataParser.set_data_directory() 
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
		WearableDataParser.set_data_directory()
		data=WearableDataParser.read_jenny_data(filename)
		
		return median(data.Steps)

	end


	res=[]
	@runAllJenny(median_steps, res)

	pl=scatter(1:length(res), res)
	display(pl)
	

end 



function getMidCompact(data, model)


	σ=2.0*π/12.0 #pm 2 hours 
	μ=HCRSimJul.ClocktoCircadianConvert(data.TimeTotal[1])
	LightProxy=HCRSimJul.interpFlux(data.TimeTotal, data.Steps)
	sol=HCRSimJul.integrateModel(LightProxy,model, [0.70, μ ∓ σ, 0.0], 0.0, 7*24.0)
	σafter= std(Array(sol)[2,end])

	return σ/σafter-1.0 


end


function genDLMOjenny(filename::String="0"; steps_multiple::Real=2.0)


	WearableDataParser.set_data_directory()
	data=WearableDataParser.read_jenny_data(subject_id=filename)
	
	 
	LightProxy=HCRSimJul.interpFlux(data.TimeTotal, data.Steps .* steps_multiple)
	model=SinglePopModel() 

	init=[0.70,HCRSimJul.ClocktoCircadianConvert(data.TimeTotal[1]), 0.0]
	tend_loop=data.TimeTotal[1]+5*24.0
	sol=integrateModel(LightProxy, model,init, data.TimeTotal[1], tend_loop)
	h=0
	while h <= 10
		sol=integrateModel(LightProxy, model, sol.u[end], data.TimeTotal[1], tend_loop)
		h+=1
	end


	DLMO=integrateObserver(LightProxy, model, sol.u[end], data.TimeTotal[1], data.TimeTotal[end], model.DLMOObs)
	if length(DLMO)>0 
		Z=1/length(DLMO)*sum(exp.( im .* DLMO .* π/12.0))
		meanDLMO=12.0/π * angle(Z) 
		if meanDLMO < 0.0 
			meanDLMO += 24.0 
		end

		compactLightScore=getMidCompact(data, model)

		idnum=data.SubjectID[1]
		@show  DLMO .% 24.0, meanDLMO
		R_DLMO=abs(Z)
		condition= (occursin("B", idnum)) ? "School" : "Summer"
		id=split(idnum, r"P|B")[1]
		return([idnum, id, condition, DLMO[end] % 24.0, compactLightScore, R_DLMO])
	else
		return missing
	end

end 



function makeJennyActogram(filename="0")

	idnum, data=WearableDataParser.read_jenny_data(filename)
	model=SinglePopModel()

	init=[0.70,HCRSimJul.ClocktoCircadianConvert(data.tstart), 0.0]
	DLMO=integrateObserver(data.Light, model, init, data.tstart, data.tend, model.DLMOObs)


	pl_actogram=HCRSimJul.makeActogram(data.Activity, data.tstart, data.tend; threshold=10.0) 
	HCRSimJul.addPhaseMarker!(pl_actogram, DLMO) 
	savefig("./ActogramsActivity/testActogramJenny_$(data.id).png")
	display(pl_actogram) 

	return 1

end

#res=[]
#@runAllJenny(makeJennyActogram,res)


function generate_file()
	res=[]
	@runAllJenny(genDLMOjenny, res)

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

	df=DataFrame(filename=fn, PatientNumber=PatientNumber, Condition=Condition, DLMO=DLMO, CompactLS=compactLightScore, DLMOPhaseAmp=R_DLMO)
	CSV.write("DLMO_Jenny_Activity_12142020_lastDLMO.csv", df)
end

# Now use this in Rprogram to fit the actual LMM


generate_file()