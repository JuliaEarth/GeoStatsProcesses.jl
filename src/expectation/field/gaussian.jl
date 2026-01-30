# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function priormean(process::GaussianProcess, domain; kwargs...)
  # process parameters
  f = process.func
  μ = process.mean
  n = nvariates(f)

  # fill domain with mean values
  vars = n > 1 ? ntuple(i -> Symbol(:field, i), n) : (:field,)
  vals = ntuple(n) do i
    fill(μ[i], nelements(domain))
  end

  # georeference mean values
  georef((; zip(vars, vals)...), domain)
end

function posteriormean(process::GaussianProcess, domain, data; kwargs...)
  # posterior mean is equivalent to simple Kriging
  model = Kriging(process.func, process.mean)
  GeoStatsModels.fitpredict(model, data, domain, neighbors=true, kwargs...)
end
