# import Pkg
# Pkg.activate(".")
include("multi_pca.jl")
using Flux, MLUtils, Random, JLD2
using Printf, Statistics, LinearAlgebra

function load_data(file)
    @load file pca_sol pca_coeff V C
    return (; V, C), (; pca_sol, pca_coeff)
end

data, pca = load_data("heat1d_pca.jld2");
println("$(size(data.V, 3)) samples loaded");

##
function convert_data(data, batchsize, train_split=0.8)
    X = reshape(data.V, :, size(data.V, 3))
    y = data.C
    (X_train, y_train), (X_test, y_test) = splitobs((X, y); at=train_split)
    return DataLoader(collect.((X_train, y_train)); batchsize, shuffle=true),
        DataLoader(collect.((X_test, y_test)); batchsize, shuffle=false)
end

train_loader, test_loader = convert_data(data, 1000);
dim_input = size(train_loader.data[1])[1]
dim_output = size(train_loader.data[2], 1)
println("input data has dimension $(dim_input)")
println("output data has dimension $(dim_output)")

##

function accuracy(loss, model, dataloader)
    tse, total = 0, 0
    for (x, y) in dataloader
        predicted = model(x)
        n = size(x, 2)
        tse += n*loss(predicted, y)
        total += n
    end
    return tse / total
end

##

model = Chain(
    Dense(dim_input => 256, relu),
    Dense(256 => 256, relu),
    Dense(256 => 256, relu),
    Dense(256 => dim_output),
)

loss = Flux.Losses.mse


# initialize the model parameters
state = Flux.setup(Flux.Adam(0.002), model)

##
for epoch in 1:20
    Flux.train!(model, train_loader, state) do m, x, y
        loss(m(x), y)
    end

    test_loss = accuracy(loss, model, test_loader)
    println("Epoch $epoch, test loss: $test_loss")
end
##

y_pred = hcat([model(x) for (x, _) in test_loader]...);
y_true = hcat([y for (_, y) in test_loader]...);
norm_diff = map(norm, eachcol(y_pred - y_true)) / 30

# fig, ax, _ = hist(log10.(norm_diff), bins=range(-3, 0.5, 31), normalization=:probability)
# ax.xlabel = "log_10(L2 error)"
# ax.ylabel = "fraction"
# fig

##
# @load "heat1d_data.jld2" diff_coeff U
coef_pred = inverse(pca.pca_coeff, y_pred)
coef_true = inverse(pca.pca_coeff, test_loader.data[2])
norm_diff = map(norm, eachcol(coef_pred - coef_true)) / 8
quantile(log10.(norm_diff), 0:.25:1)
