"""MIT License

Copyright (c) [2025] [Nikki Walker]

File that uses MC class"""
import montecarlo_class as MC


# Test Case 1: n-dimension sphere

# Initialise MC in 2D
MC_rand_2D = MC.MC(d=2, n=100000000, seed=134)

# Calculate exact volume for D case
volume_2D_exact = MC_rand_2D.volume_exact()

# Calculate estimated volume using MC
volume_2D_MC = MC_rand_2D.volume(r=1)
volume_2D_MC_mean, volume_2D_MC_err, volume_2D_variance = MC_rand_2D.statistics(volume_2D_MC)


if MC_rand_2D.rank == 0:
    print("Test Case 1: n-dimension sphere", "\n")
    print("2D sphere:", "\n")
    print(f"Exact volume: {volume_2D_exact:.6f}")
    print(f"MC estimate: {volume_2D_MC:.6f} +/- {volume_2D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_2D_MC - volume_2D_exact)/volume_2D_exact*100:.2f}%")
    MC_rand_2D.timings()

# Initialise MC in 3D
MC_rand_3D = MC.MC(d=3, n=100000000, seed=134)

# Calculate exact volume for D case
volume_3D_exact = MC_rand_3D.volume_exact()

# Calculate estimated volume using MC
volume_3D_MC = MC_rand_3D.volume(r=1)
volume_3D_MC_mean, volume_3D_MC_err, volume_3D_variance = MC_rand_3D.statistics(volume_3D_MC)

if MC_rand_3D.rank == 0:
    print("3D sphere:", "\n")
    print(f"Exact volume: {volume_3D_exact:.6f}")
    print(f"MC estimate: {volume_3D_MC:.6f} +/- {volume_3D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_3D_MC - volume_3D_exact)/volume_3D_exact*100:.2f}%")
    MC_rand_3D.timings()

# Initialise MC in 4D
MC_rand_4D = MC.MC(d=4, n=100000000, seed=134)

# Calculate exact volume _for D case
volume_4D_exact = MC_rand_4D.volume_exact()

# Calculate estimated volume using MC
volume_4D_MC = MC_rand_4D.volume(r=1)
volume_4D_MC_mean, volume_4D_MC_err, volume_4D_variance = MC_rand_4D.statistics(volume_4D_MC)

if MC_rand_4D.rank == 0:
    print("4D sphere:", "\n")
    print(f"Exact volume: {volume_4D_exact:.6f}")
    print(f"MC estimate: {volume_4D_MC:.6f} +/- {volume_4D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_4D_MC - volume_4D_exact)/volume_4D_exact*100:.2f}%")
    MC_rand_4D.timings()

# Initialise MC in 5D
MC_rand_5D = MC.MC(d=5, n=100000000, seed=134)

# Calculate exact volume for D case
volume_5D_exact = MC_rand_5D.volume_exact()

# Calculate estimated volume using MC
volume_5D_MC = MC_rand_5D.volume(r=1)
volume_5D_MC_mean, volume_5D_MC_err, volume_5D_variance = MC_rand_5D.statistics(volume_5D_MC)

if MC_rand_5D.rank == 0:
    print("5D sphere:", "\n")
    print(f"Exact volume: {volume_5D_exact:.6f}")
    print(f"MC estimate: {volume_5D_MC:.6f} +/- {volume_5D_MC_err:.10f}")
    print(f"Relative error: {abs(volume_5D_MC - volume_5D_exact)/volume_5D_exact*100:.2f}%")
    MC_rand_5D.timings()

# Test Case 2: Evaluation of gaussian function

# Initialise MC for 1D
MC_rand_1D = MC.MC(d=1, n=100000000, seed=134)

#Calculate MC estimate for 1D gaussian for range of values
MC_1D_gauss_val1 = MC_rand_1D.gauss_integ(sigma=[1.0], x0=[0])

if MC_rand_1D.rank == 0:
    print("Test Case 2: Gaussian Integration: \n")
    print("1D Gaussian sigma = 1.0 and x0 = 0", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_1D_gauss_val1[0]:.6f} +/- {MC_1D_gauss_val1[1]:.6f}")
    print(f"Variance: {MC_1D_gauss_val1[2]:.6f}")
    MC_rand_1D.timings()


MC_1D_gauss_val2 = MC_rand_1D.gauss_integ(sigma=[2.0], x0=[2.0])

if MC_rand_1D.rank == 0:
    print("1D Gaussian sigma = 2.0 and x0 = 2.0", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_1D_gauss_val2[0]:.6f} +/- {MC_1D_gauss_val2[1]:.6f}")
    print(f"Variance: {MC_1D_gauss_val2[2]:.6f}")
    MC_rand_1D.timings()


MC_1D_gauss_val3 = MC_rand_1D.gauss_integ(sigma=[3.0], x0=[1.0])

if MC_rand_1D.rank == 0:
    print("1D Gaussian sigma = 3.0 and x0 = 1.0", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_1D_gauss_val3[0]:.6f} +/- {MC_1D_gauss_val3[1]:.6f}")
    print(f"Variance: {MC_1D_gauss_val3[2]:.6f}")
    MC_rand_1D.timings()

# Initialise MC for 6D
MC_rand_6D = MC.MC(d=6, n=100000000, seed=134)

#Calculate MC estimate for 6D gaussian
MC_6D_gauss_val1 = MC_rand_6D.gauss_integ(sigma=[1.0, 1.0, 2.0, 1.0, 2.0, 2.0],
x0=[0, 1, 0, 1, 0, 1])

if MC_rand_6D.rank == 0:
    print("6D Gaussian sigma = [1.0, 1.0, 2.0, 1.0, 2.0, 2.0] and x0 = [0, 1, 0, 1, 0, 1]", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_6D_gauss_val1[0]:.6f} +/- {MC_6D_gauss_val1[1]:.6f}")
    print(f"Variance: {MC_6D_gauss_val1[2]:.6f}")
    MC_rand_6D.timings()


MC_6D_gauss_val2 = MC_rand_6D.gauss_integ(sigma=[3.0, 2.0, 3.0, 1.0, 3.0, 1.0],
x0=[-1, 0, 1, -1, -1, 1])

if MC_rand_6D.rank == 0:
    print("6D Gaussian sigma = [3.0, 2.0, 3.0, 1.0, 3.0, 1.0] and x0 = [-1, 0, 1, -1, -1, 1]", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_6D_gauss_val2[0]:.6f} +/- {MC_6D_gauss_val2[1]:.6f}")
    print(f"Variance: {MC_6D_gauss_val2[2]:.6f}")
    MC_rand_6D.timings()


MC_6D_gauss_val3 = MC_rand_6D.gauss_integ(sigma=[2.0, 3.0, 3.0, 3.0, 1.0, 2.0],
x0=[1, -1, 0, 1, 0, 0])

if MC_rand_6D.rank == 0:
    print("6D Gaussian sigma = [2.0, 3.0, 3.0, 3.0, 1.0, 2.0] and x0 = [1, -1, 0, 1, 0, 0]", "\n")
    print("Exact Gaussian: 1.0")
    print(f"MC estimate: {MC_6D_gauss_val3[0]:.6f} +/- {MC_6D_gauss_val3[1]:.6f}")
    print(f"Variance: {MC_6D_gauss_val3[2]:.6f}")
    MC_rand_6D.timings()
