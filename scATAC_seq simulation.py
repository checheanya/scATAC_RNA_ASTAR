import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import nbinom

def simulate_single_cell_atac_seq(num_cells, num_regions, library_size, mean_read_count, dispersion):
    # Generate the library size for each cell (drawn from a normal distribution)
    cell_library_sizes = np.random.normal(loc=library_size, scale=library_size/2, size=num_cells).astype(int)
    
    # Simulate read counts for each cell and region
    simulated_data = np.zeros((num_cells, num_regions), dtype=int)
    for cell_idx in range(num_cells):
        # Simulate read counts for each region in the cell
      
        for region_idx in range(num_regions):
          
            # NEGATIVE BINOMIAL
            read_count = nbinom.rvs(n=dispersion, p=dispersion / (dispersion + mean_read_count))
            # POISSON
            # read_count = poisson.rvs(mean_read_count)
            # BERNOULLI
            # read_count = bernoulli.rvs(mean_read_count / library_size)
            # BINOMIAL
            # read_count = binom.rvs(library_size, mean_read_count / library_size)
            # ZERO-INFLATED NB
            '''
            zero_prob = 0.2  # Probability of a zero count
                  zero_count = np.random.choice([0, 1], p=[1 - zero_prob, zero_prob])
                  if zero_count == 0:
                      read_count = nbinom.rvs(n=dispersion, p=dispersion / (dispersion + mean_read_count))
                  else:
                      read_count = 0
            '''

            # Ensure the read count does not exceed the library size
            read_count = min(read_count, cell_library_sizes[cell_idx])
            # Assign the read count to the corresponding cell and region
            simulated_data[cell_idx, region_idx] = read_count
            # Subtract the assigned read count from the library size of the cell
            cell_library_sizes[cell_idx] -= read_count
    
    return simulated_data

# Parameters
num_cells = 100
num_regions = 5000
library_size = 10000
mean_read_count = 5
dispersion = 1

# Simulate data
simulated_data = simulate_single_cell_atac_seq(num_cells, num_regions, library_size, mean_read_count, dispersion)

# Generate plots to illustrate the distributions
plt.figure(figsize=(15, 6))

# Plot the Negative Binomial distribution
plt.subplot(1, 5, 1)
nbinom_data = nbinom.rvs(n=dispersion, p=dispersion / (dispersion + mean_read_count), size=num_regions)
plt.hist(nbinom_data, bins=30, density=True, color='blue', alpha=0.6)
plt.title('Negative Binomial')
plt.xlabel('Read Count')
plt.ylabel('Density')

# Plot the Poisson distribution
plt.subplot(1, 5, 2)
poisson_data = poisson.rvs(mean_read_count, size=num_regions)
plt.hist(poisson_data, bins=30, density=True, color='green', alpha=0.6)
plt.title('Poisson')
plt.xlabel('Read Count')
plt.ylabel('Density')

# Plot the Zero-Inflated Negative Binomial distribution
plt.subplot(1, 5, 3)
zero_prob = 0.2
zinf_data = []
for _ in range(num_regions):
    zero_count = np.random.choice([0, 1], p=[1 - zero_prob, zero_prob])
    if zero_count == 0:
        zinf_data.append(nbinom.rvs(n=dispersion, p=dispersion / (dispersion + mean_read_count)))
    else:
        zinf_data.append(0)
plt.hist(zinf_data, bins=30, density=True, color='red', alpha=0.6)
plt.title('Zero-Inflated Negative Binomial')
plt.xlabel('Read Count')
plt.ylabel('Density')

# Plot the Bernoulli distribution
plt.subplot(1, 5, 4)
bernoulli_data = bernoulli.rvs(mean_read_count / library_size, size=num_regions)
plt.hist(bernoulli_data, bins=30, density=True, color='green', alpha=0.6)
plt.title('Bernoulli')
plt.xlabel('Read Count')
plt.ylabel('Density')

# Plot the Binomial distribution
plt.subplot(1, 5, 5)
binom_data = binom.rvs(library_size, mean_read_count / library_size, size=num_regions)
plt.hist(binom_data, bins=30, density=True, color='purple', alpha=0.6)
plt.title('Binomial')
plt.xlabel('Read Count')
plt.ylabel('Density')


plt.tight_layout()
plt.show()
