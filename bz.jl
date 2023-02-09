using DSP, Plots

dim = 200
frames = 150
diffusion = 5
reaction_rates = [[1.0, 1.0],
                  [1.0, 1.0],
                  [1.0, 1.0]]

mesh = [rand(Float64, (dim, dim)), rand(Float64, (dim, dim)), rand(Float64, (dim, dim))]
dmesh = similar(mesh)
gaussian_distribution = ℯ.^(-(range(-ℯ, ℯ, diffusion*2+1).^2))
diffusion_kernel = gaussian_distribution/sum(gaussian_distribution)

function diffuse(reactant)
    return conv(diffusion_kernel, diffusion_kernel, reactant)[1+diffusion:end-diffusion, 1+diffusion:end-diffusion]
end

@gif for i=1:frames
    global mesh, demesh
    mesh .= diffuse.(mesh)
    for i=1:3  dmesh[i] = mesh[i] .* ( reaction_rates[i][1] * mesh[(i)%3+1] .- reaction_rates[i][2] * mesh[(i+1)%3+1] )  end
    mesh .+= dmesh
    clamp!.(mesh, 0, 1)
    heatmap(mesh[1],  size=(440, 400))
end
