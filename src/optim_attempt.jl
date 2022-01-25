using Revise
using tbconvax, BenchmarkTools, CSV, Tables, Optim, YAML, Dates
function f(x)
    try
        sobgen(target_extract(main(fixed_input = fi, variable_input = pst(x))))
    catch
        return 300
    end
end

CMA = YAML.load_file("data/CMA_indicator.yml")["CMA"]
lower = zeros(Float64, 25);
upper = ones(Float64, 25);
initial_x = rand(25);

opt = optimize(
    f,
    lower,
    upper,
    initial_x,
    Fminbox(NelderMead()),
    Optim.Options(show_trace = true, show_every = 500),
)

stamp = Dates.format(Dates.now(), Dates.DateFormat("yyyy-mm-dd-HHMM"))

id = "$(stamp)_$(CMA)_seed"

YAML.write_file(
    joinpath("output", CMA, "params", "seeds", "$stamp.yml"),
    Dict("stamp" => stamp, "id" => id, "ts" => ts),
)
YAML.write_file(
    joinpath("output", CMA, "params", "seeds", "latest.yml"),
    Dict("stamp" => stamp, "id" => id, "ts" => ts),
)
