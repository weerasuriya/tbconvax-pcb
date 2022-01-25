using Revise
using tbconvax, BenchmarkTools, CSV, Tables, Optim, YAML, Dates
inc2050hm = YAML.load_file("output/HM/targets/inc2050/latest.yml")

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

ts = [
    0.725393269552202,
    0.885535185349462,
    0.6596469781229015,
    0.8204147664115499,
    3.3999114416094374e-5,
    0.9355613628569176,
    0.47274157421109636,
    0.9582149797660238,
    0.9387010381708407,
    0.9371828094744447,
    0.0459094279066223,
    0.9859142436024383,
    0.15260223942134882,
    0.02394575944272129,
    0.16242098959646323,
    0.6763093472778751,
    0.2522324576553068,
    0.5889543364697994,
    2.87841532496171e-5,
    0.4480121705626443,
    1.5340864349697528e-5,
    0.2932057994629847,
    0.02356621573251142,
    0.07971094176691093,
    0.2245433161650628,
]

id = "$(stamp)_$(CMA)_seed"

YAML.write_file(
    joinpath("output", CMA, "params", "seeds", "$stamp.yml"),
    Dict("stamp" => stamp, "id" => id, "ts" => ts, "inc2050hm_tgt" => inc2050hm),
)
YAML.write_file(
    joinpath("output", CMA, "params", "seeds", "latest.yml"),
    Dict("stamp" => stamp, "id" => id, "ts" => ts, "inc2050hm_tgt" => inc2050hm),
)
