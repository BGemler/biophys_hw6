import matplotlib.pyplot as plt
import random
import csv
import math
from numpy import dot
import numpy as np

random.seed(2021) # fix for reproducibility 

def perform_polymer_walk(N, alpha):
	"""
	note: python default conducts sin/cos in RADIANS
	"""
	positions = [[0, 0]] # [x, y]
	phi = 0

	for i in range(1, N - 1):
		# pick new angle from -pi to pi, scaled by alpha
		deltaphi = alpha * random.uniform(-1, 1) * math.pi

		phi = phi + deltaphi

		deltax = math.cos(phi)
		deltay = math.sin(phi)

		# set x, y positions at step
		new_x = positions[i-1][0] + deltax
		new_y = positions[i-1][1] + deltay

		positions.append([new_x, new_y])

	return positions


def compute_persist_length(positions, max_k, N):
	"""
	persistance length calculated for each k
	"""
	# first, need to compute the segment vectors
	segment_vectors = []
	for i in range(N-2):
		x_2, y_2 = positions[i + 1]
		x_1, y_1 = positions[i]

		segment_vectors.append([x_2 - x_1, y_2 - y_1])

	# compute c(k) for k in range 1, ... 199)
	all_cks = []
	for k in range(1, max_k):
		# first, compute the summation term
		sum_term = 0.0
		for i in range(N - max_k - 1):
			segment_vector_i = np.asarray(segment_vectors[i])
			segment_vector_i_k = np.asarray(segment_vectors[i + k])

			# numpy dot product 
			dot_product = segment_vector_i.dot(segment_vector_i_k)

			sum_term += dot_product

		# find c(k) at this k value
		c_k = (1 / (N - max_k)) * sum_term

		all_cks.append([k, c_k])

	return all_cks


def plot_polymer_walk(positions, N, alpha):
	"""
	"""
	out_img = "results/polymer-walk-alpha=" + str(alpha) + ".png"

	plot_x, plot_y = [], []
	for x, y in positions:
		plot_x.append(x)
		plot_y.append(y)

	plt.scatter(plot_x, plot_y)
	plt.xlabel("X Position")
	plt.ylabel("Y Position")

	plt.title("Polymer Walk for Alpha=" + str(alpha) + "\n" + \
								"Number of Steps (N)=" + str(N))

	plt.savefig(out_img, bbox_inches='tight')
	plt.close()

	return


def plot_ck_values(all_cks, alpha, N, max_k, log_plot):
	"""
	"""
	out_img = "results/ck-values-log=" + str(log_plot) + \
								"-alpha=" + str(alpha) + ".png"

	plot_x, plot_y = [], []
	for x, y in all_cks:
		plot_x.append(x)
		plot_y.append(y)

	#if log_plot == False:
	plt.scatter(plot_x, plot_y)
	if log_plot == True:
		plt.yscale('log')

	plt.xlabel("k")
	plt.ylabel("c(k)")

	plt.title("c(k) for k in range (1," + str(max_k) + ")" + "\n" + \
								"Y-Axis Log?" + str(log_plot) + "\n" + \
								"Alpha=" + str(alpha) + "\n" + \
								"Number of Steps (N)=" + str(N))

	plt.savefig(out_img, bbox_inches='tight')
	plt.close()

	if log_plot == False:
		with open("results/ck-values.csv", "w") as f:
			out = csv.writer(f)
			out.writerow(["X", "Y"])

			for x, y in all_cks:
				out.writerow([x, y])
		f.close()


	return


def main():
	"""
	perform the polymer plot, adjusting variables listed below
	"""
	N = 10000
	alpha = 0.05
	max_k = 200

	# get polymer walk positions 
	positions = perform_polymer_walk(N, alpha)

	# plot the polymer walk
	plot_polymer_walk(positions, N, alpha)

	# get c(k) values
	all_cks = compute_persist_length(positions, max_k, N)

	# plot c(k) values
	plot_ck_values(all_cks, alpha, N, max_k, log_plot = False)
	plot_ck_values(all_cks, alpha, N, max_k, log_plot = True)

	return


main()