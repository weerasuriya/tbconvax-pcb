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
    0.9055246373740291,
    0.9105295412943408,
    0.9529431259760667,
    0.9903608779207261,
    0.20986647119791696,
    0.09944145907330644,
    0.4286313099528051,
    0.8644011392670478,
    0.7525016735770064,
    0.24644498560916409,
    0.019291215791808404,
    0.8156941077351216,
    0.02868053933730207,
    0.020944142449265957,
    0.047889671755414545,
    0.5745265941373922,
    0.9416406642657684,
    0.9149341908848218,
    0.41368836466631353,
    0.8294255451744911,
    0.8821607330946964,
    0.9193243553026135,
    0.10198799405152513,
    0.6929062236740211,
    0.32471459235791283,
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
