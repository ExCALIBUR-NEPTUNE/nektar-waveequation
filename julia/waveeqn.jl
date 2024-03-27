using FFTW, LinearAlgebra, Plots

function waveeqn()
  N = 64 # number of cells in x and y
  L = 1.0 # length in x and y
  x = y = collect(1/2N:1/N:1-1/2N) .* L

  v = 1 # the advection speed
  dx = L / N # linear cell size
  dt = dx / v * 2# * 0.45 # timestep
  k0 = 2pi/L # small wavenumber
  function fxy(x, y, t) # initial conditions
    kx = ky = k0
    ω = sqrt(kx^2 + ky^2) * v
    return sin(kx * x + ky * y - ω * t)
  end

  # hold the fields in these at n+1 (z1), n (z0) and n-1 (z_1)
  z_1 = fft([fxy(xi, yi, -dt) for xi in x, yi in y])
  z0 = fft([fxy(xi, yi, 0.0) for xi in x, yi in y])
  z1 = zeros(ComplexF64, size(z0)...)

  # 1d wavenumbers in the order of the FFT
  k1d = im .* [0:N÷2-1; -N÷2:-1] .* k0
  # 2d waveumber squared matrix
  ∇² = [kx^2 + ky^2 for kx in k1d, ky in k1d]
  θ = 0.5 # the implicitness parameter (0.5 for unconditionally implicit)

  NT = 2L÷(v * dt) # number of steps to solve for
  @info "Running for $NT timesteps"
  anim = @animate for tstep in 1:NT
    if iszero(θ) # explicit
      @. z1 = (2 + dt^2 * ∇²) * z0 - z_1
    else # implicit
      λ = 2 / θ / dt^2
      @. z1 = (-2 * (λ + (1-θ)/θ * ∇²) * z0 - (∇² - λ) * z_1) / (∇² - λ)
    end
    z_1 .= z0 # copy z(n) to z(n-1)
    z0 .= z1 # copy z(n+1) to z(n)
    heatmap(real.(ifft(z0))) # plot for the gif
  end
  gif(anim, "waveeqn.jl.gif", fps=16)
end

waveeqn() # run this thing


