# import Pkg
# Pkg.activate(".")
using Flux, MLUtils, Random, JLD2
using Printf, Statistics, LinearAlgebra, MultivariateStats

function load_data(file)
    @load file Coeff_k Coeff_c U
    return (; U, K=Coeff_k, C=Coeff_c)
end

data = load_data("advdiff1d_bc_data2.jld2");
println("$(size(data.U, 3)) samples loaded");

##
function convert_data(data, batchsize, train_split=0.8)
    X = reshape(data.U, size(data.U, 1), :)
    u_pod = fit(PCA, X; maxoutdim=30, pratio=0.99999)
    V = reshape(predict(u_pod, X), :, size(data.U)[2:3]...)
    X = reshape(permutedims(V, (2, 1, 3)), size(V, 2), :)
    u_pca = fit(PCA, X; maxoutdim=30, pratio=0.99999)
    X = reshape(predict(u_pca, X), :, size(V,3))
    c_pca = fit(PCA, data.C; maxoutdim=11, pratio=0.99999)
    c = predict(c_pca, data.C)
    k_pca = fit(PCA, data.K; maxoutdim=11, pratio=0.99999)
    k = predict(k_pca, data.K)
    y = vcat(k, c)
    (X_train, y_train), (X_test, y_test) = splitobs((X, y); at=train_split)
    return (; ux=u_pod, ut=u_pca, c=c_pca, k=k_pca),
     DataLoader(collect.((X_train, y_train)); batchsize, shuffle=true),
        DataLoader(collect.((X_test, y_test)); batchsize, shuffle=false)
end

pca, train_loader, test_loader = convert_data(data, 200);
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

model = Chain(
    Dense(dim_input => 400, tanh),
    Dense(400 => 200, tanh),
    Dense(200 => dim_output),
)

loss = Flux.Losses.mse
state = Flux.setup(Flux.Adam(0.001), model)

##
for epoch in 1:200
    Flux.train!(model, train_loader, state) do m, x, y
        loss(m(x), y)
    end

    test_loss = accuracy(loss, model, test_loader)
    println("Epoch $epoch, test loss: $test_loss")
end
##

function coeffs(y, pca)
    n = size(pca.k.prinvars, 1)
    k = reconstruct(pca.k, y[1:n, :])
    c = reconstruct(pca.c, y[n+1:end, :])
    return k, c
end

y_pred = hcat([model(x) for (x, _) in test_loader]...);
k_pred, c_pred = coeffs(y_pred, pca)
y_true = hcat([y for (_, y) in test_loader]...);
k_true, c_true = coeffs(test_loader.data[2], pca)
norm_diff_c = map(norm, eachcol(c_pred - c_true)) / 8
norm_diff_k = map(norm, eachcol(k_pred - k_true)) / 8
qc = quantile(log10.(norm_diff_c), 0:.25:1)
qk = quantile(log10.(norm_diff_k), 0:.25:1);

##
fig = Figure(colormap="matter", markersize=5)
ntest = size(y_true, 2)
axkk = Axis(fig[1, 1], xlabel="min k", ylabel="error in k")
scatter!(map(minimum, eachcol(k_true)), norm_diff_k, color=1:ntest)
axkc = Axis(fig[2, 1], xlabel="min k", ylabel="error in c")
scatter!(map(minimum, eachcol(k_true)), norm_diff_c, color=1:ntest)
axck = Axis(fig[1, 2], xlabel="min c", ylabel="error in k")
scatter!(map(minimum, eachcol(c_true)), norm_diff_k, color=1:ntest)
axcc = Axis(fig[2, 2], xlabel="min c", ylabel="error in c")
scatter!(map(minimum, eachcol(c_true)), norm_diff_c, color=1:ntest)
fig
