# import Pkg
# Pkg.activate(".")
using Flux, MLUtils, Random, JLD2
using Printf, Statistics, LinearAlgebra, MultivariateStats

function load_data(file)
    @load file Coeff_kx Coeff_ky U
    # return (; U=U[:,:,:,1:950], Kx=Coeff_kx[:,:,1:950], Ky=Coeff_ky[:,:,1:950])
    return (; U=U, Kx=Coeff_kx, Ky=Coeff_ky)
end

data = load_data("heat2d_data3.jld2");
@load "heat2d_data3.jld2" x t
nx = length(x)
println("$(size(data.U, 4)) samples loaded");

##
function convert_data(data, batchsize, train_split=0.8)
    X = reshape(data.U, prod(size(data.U)[1:2]), :)
    u_pod = fit(PCA, X; maxoutdim=200, pratio=0.9999)
    V = reshape(predict(u_pod, X), :, size(data.U)[3:4]...)
    X = reshape(permutedims(V, (2, 1, 3)), size(V, 2), :)
    u_pca = fit(PCA, X; maxoutdim=30, pratio=0.99999)
    X = reshape(predict(u_pca, X), :, size(V,3))

    n = size(data.Kx, 3)
    K = reshape(data.Kx, :, n)
    kx_pca = fit(PCA, K; maxoutdim=40, pratio=0.99999)
    kx = predict(kx_pca, K)
    K = reshape(data.Ky, :, n)
    ky_pca = fit(PCA, K; maxoutdim=40, pratio=0.99999)
    ky = predict(ky_pca, K)
    y = vcat(kx, ky)
    (X_train, y_train), (X_test, y_test) = splitobs((X, y); at=train_split)
    return (; uxy=u_pod, ut=u_pca, kx=kx_pca, ky=ky_pca),
     DataLoader(collect.((X_train, y_train)); batchsize, shuffle=true),
        DataLoader(collect.((X_test, y_test)); batchsize, shuffle=false)
end

pca, train_loader, test_loader = convert_data(data, 100);
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
    Dense(dim_input => 600, tanh),
    Dense(600 => 400, tanh),
    Dense(400 => 200, tanh),
    Dense(200 => dim_output),
)

loss = Flux.Losses.mse
state = Flux.setup(Flux.Adam(0.0005), model)

##
for epoch in 1:100
    Flux.train!(model, train_loader, state) do m, x, y
        loss(m(x), y)
    end

    test_loss = accuracy(loss, model, test_loader)
    println("Epoch $epoch, test loss: $test_loss")
end
##

function coeffs(y, pca)
    N = size(pca.kx.prinvars, 1)
    k = reconstruct(pca.kx, y[1:N, :])
    Kx = reshape(k, nx, nx, :)
    k = reconstruct(pca.ky, y[N+1:end, :])
    Ky = reshape(k, nx, nx, :)
    return Kx, Ky
end

y_pred = hcat([model(x) for (x, _) in test_loader]...);
Kx_pred, Ky_pred = coeffs(y_pred, pca)
y_true = hcat([y for (_, y) in test_loader]...);
Kx_true, Ky_true = coeffs(test_loader.data[2], pca)

scale = hypot.(map(norm, eachslice(Kx_true, dims=3)), map(norm, eachslice(Ky_true, dims=3)))
diff_kx = map(norm, eachslice(Kx_pred - Kx_true, dims=3)) ./ scale
diff_ky = map(norm, eachslice(Ky_pred - Ky_true, dims=3)) ./ scale
(quantile(log10.(diff_kx), 0:.25:1), quantile(log10.(diff_ky), 0:.25:1))
hist([diff_kx;diff_ky])

##
# fig = Figure(colormap="matter", markersize=5)
# ntest = size(y_true, 2)
# axkk = Axis(fig[1, 1], xlabel="min k", ylabel="error in k")
# scatter!(map(minimum, eachcol(k_true)), norm_diff_k, color=1:ntest)
# axkc = Axis(fig[2, 1], xlabel="min k", ylabel="error in c")
# scatter!(map(minimum, eachcol(k_true)), norm_diff_c, color=1:ntest)
# axck = Axis(fig[1, 2], xlabel="min c", ylabel="error in k")
# scatter!(map(minimum, eachcol(c_true)), norm_diff_k, color=1:ntest)
# axcc = Axis(fig[2, 2], xlabel="min c", ylabel="error in c")
# scatter!(map(minimum, eachcol(c_true)), norm_diff_c, color=1:ntest)
# fig
