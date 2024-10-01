# import Pkg
# Pkg.activate(".")
using Flux, MLUtils, Random, JLD2
using Printf, Statistics, LinearAlgebra, MultivariateStats

function load_data(file)
    @load file diff_coeff U
    return (; U, coeff=diff_coeff)
end

data = load_data("heat1d_data.jld2");
println("$(size(data.U, 3)) samples loaded");

##
function convert_data(data, batchsize, train_split=0.8)
    X = reshape(data.U, size(data.U, 1), :)
    u_pod = fit(PCA, X; maxoutdim=30, pratio=0.9999)
    V = reshape(predict(u_pod, X), :, size(data.U)[2:3]...)
    X = reshape(permutedims(V, (2, 1, 3)), size(V, 2), :)
    u_pca = fit(PCA, X; maxoutdim=30, pratio=0.9999)
    X = reshape(predict(u_pca, X), :, size(V,3))
    c_pca = fit(PCA, data.coeff; maxoutdim=15, pratio=0.9999)
    y = predict(c_pca, data.coeff)
    (X_train, y_train), (X_test, y_test) = splitobs((X, y); at=train_split)
    return (; u_pod, u_pca, c_pca),
     DataLoader(collect.((X_train, y_train)); batchsize, shuffle=true),
        DataLoader(collect.((X_test, y_test)); batchsize, shuffle=false)
end

pca, train_loader, test_loader = convert_data(data, 1000);
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
    Dense(dim_input => 128, tanh),
    Dense(128 => 128, tanh),
    Dense(128 => dim_output),
)

loss = Flux.Losses.mse
# initialize the model parameters
state = Flux.setup(Flux.Adam(0.001), model)

##
for epoch in 1:200
    Flux.train!(model, train_loader, state) do m, x, y
        loss(m(x), y)
    end

    # test_loss = accuracy(loss, model, test_loader)
    # println("Epoch $epoch, test loss: $test_loss")
end
##

y_pred = hcat([model(x) for (x, _) in test_loader]...);
y_true = hcat([y for (_, y) in test_loader]...);
coef_pred = reconstruct(pca.c_pca, y_pred)
coef_true = reconstruct(pca.c_pca, test_loader.data[2])
norm_diff = map(norm, eachcol(coef_pred - coef_true)) / 8
qq = quantile(log10.(norm_diff), 0:.25:1)
