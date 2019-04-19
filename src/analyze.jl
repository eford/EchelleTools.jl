using Statistics

function mean_spectrum(data::AbstractEchelleObservationSet)
   sum(data.flux./data.var,dims=3)./sum(one(eltype(data.var))./data.var,dims=3)
end

function median_spectrum(data::AbstractEchelleObservationSet)
   median(data.flux,dims=3)
end


