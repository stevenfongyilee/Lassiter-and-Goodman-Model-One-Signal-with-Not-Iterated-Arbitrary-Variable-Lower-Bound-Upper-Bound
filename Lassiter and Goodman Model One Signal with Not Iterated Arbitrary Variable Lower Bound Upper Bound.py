# Here we allow arbitrary iteration up the hierarchy.

import numpy
numpy.set_printoptions(linewidth = 120)
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import uniform
from scipy.stats import truncnorm
from scipy import integrate

###################

def receiver_0_signal_1(h, theta):
	if h < theta:
		return 0.
	else:
		return state_distribution.pdf(h) / (1. - state_distribution.cdf(theta))
		
def receiver_0_signal_not_1(h, theta):
	if h > theta:
		return 0.
	else:
		return state_distribution.pdf(h) / (state_distribution.cdf(theta))

def pragmatic_sender(n):
	if n not in pragmatic_sender_memo:
		pragmatic_receiver_signal_1_theta_1_by_h_array = pragmatic_receiver(n - 1)[0]
		pragmatic_receiver_signal_not_1_theta_1_by_h_array = pragmatic_receiver(n - 1)[1]
		
# note that we have two options here: first we could make the pragmatic sender sensitive
# to just h, second we could make it sensitive to h and theta. of course in Lassiter and
# Goodman's model the first level sender is only sensitive to h. if we want to keep a
# single model all the way up the hierarchy it seems like making the sender only h
# sensitive is the way to go. there is also a question about theta values: do we get them
# from the same level or the previous level? What way to think about it is consistent with
# the initial model? does the initial model (the first level sender) prescribe only one
# way of doing the calculations? Perhaps one way is to think of the receiver at any level
# as taking the background knowledge of h and then applying a literal interpretation which
# simply renormalizes above any given theta. at the bottom level 0 the background
# information is an absolute prior, and at level 2 the background for any given signal is
# the h posterior. Then the level 3 sender thinks of the level 2 receiver as renormalizing
# for any given theta above (or below) the h posterior for the given signal. That can give
# rise to a level 4 receiver who in turn derives a new posterior. Is this any different
# than what we have below for h sensitive?
		
		if pragmatic_sender_type == 'h_sensitive_version_3':
		
			pragmatic_receiver_signal_1_h_array = numpy.sum(pragmatic_receiver(n - 1)[0], axis = (0))
			pragmatic_receiver_signal_not_1_h_array = numpy.sum(pragmatic_receiver(n - 1)[1], axis = (0))
						
			pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.tile(pragmatic_receiver_signal_1_h_array / (num_states), (num_states, 1))) - cost))
			pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.tile(pragmatic_receiver_signal_not_1_h_array / (num_states), (num_states, 1))) - (cost + cost_of_not)))
						
			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
			
			denominator_array = pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized
		
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		elif pragmatic_sender_type == 'h_sensitive':
		
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[0], axis = (1)), (num_states, 1)).T
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[1], axis = (1)), (num_states, 1)).T
			
			print 'pragmatic_receiver_signal_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_1_theta_1_by_h_array
			print 'pragmatic_receiver_signal_not_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_by_h_array
						
			print 'pragmatic_receiver_signal_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_1_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_not_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_not_1_pre_normalization_factor_array
						
			pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array) - cost))
			pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array) - (cost + cost_of_not)))

			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		elif pragmatic_sender_type == 'h_and_theta_sensitive':
		
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.ones((num_states, num_states))
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.ones((num_states, num_states))
			
			pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array) - cost))
			pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array) - (cost + cost_of_not)))

			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)
		
		elif pragmatic_sender_type == 'modified h sensitive':
		
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[0], axis = (1)), (num_states, 1)).T
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[1], axis = (1)), (num_states, 1)).T
			
			print 'pragmatic_receiver_signal_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_1_theta_1_by_h_array
			print 'pragmatic_receiver_signal_not_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_by_h_array
						
			print 'pragmatic_receiver_signal_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_1_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_not_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_not_1_pre_normalization_factor_array
						
			pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.sum(weighting_array.reshape(1, num_states, num_states) * (pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array).reshape(num_states, 1, num_states), axis = 2)) - cost))
			pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.sum(weighting_array.reshape(1, num_states, num_states) * (pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array).reshape(num_states, 1, num_states), axis = 2)) - (cost + cost_of_not)))

			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		pragmatic_sender_memo[n] = numpy.asarray((pragmatic_sender_signal_1_normalized, pragmatic_sender_signal_not_1_normalized))
	return pragmatic_sender_memo[n]

def pragmatic_receiver(n):
	if n not in pragmatic_receiver_memo:

		if theta_posterior_source == 'signal 1':
			theta_1_distribution_array = numpy.sum(pragmatic_receiver(n-2)[0:1], axis = (0, 2))
		elif theta_posterior_source == 'signal 1, not 1':
			theta_1_distribution_array = (numpy.sum(pragmatic_receiver(n-2)[0:], axis = (0, 2)) / 2.)
		elif theta_posterior_source == 'signal_specific':
			theta_1_distribution_array_array = numpy.sum(pragmatic_receiver(n-2), axis = 2)
			
		if theta_posterior_source != 'signal_specific':
			theta_1_distribution_array_array = numpy.tile(theta_1_distribution_array, (2, 1))
		
		for array_num in range(len(theta_1_distribution_array_array)):

			plot_theta_1_distribution_array = theta_1_distribution_array_array[array_num]
			
			fig, ax = pyplot.subplots(1,1)
			pyplot.plot(array_0, plot_theta_1_distribution_array)
			pyplot.show()
			pyplot.close()

		theta_1_distribution_array_array = numpy.swapaxes(theta_1_distribution_array_array, 0, 1)
		
		pragmatic_receiver_signal_1_h_prior = numpy.sum(pragmatic_receiver(n-2)[0], axis = (0))
		pragmatic_receiver_signal_not_1_h_prior = numpy.sum(pragmatic_receiver(n-2)[1], axis = (0))
		
		if n <= 2:
			pragmatic_receiver_signal_1_pre_normalized = pragmatic_sender(n-1)[0] * theta_1_distribution_array_array[:,0:1] * pragmatic_receiver_signal_0_h_prior
			pragmatic_receiver_signal_not_1_pre_normalized = pragmatic_sender(n-1)[1] * theta_1_distribution_array_array[:,1:2] * pragmatic_receiver_signal_0_h_prior
		elif n > 2:
			pragmatic_receiver_signal_1_pre_normalized = pragmatic_sender(n-1)[0] * theta_1_distribution_array_array[:,0:1] * pragmatic_receiver_signal_1_h_prior
			pragmatic_receiver_signal_not_1_pre_normalized = pragmatic_sender(n-1)[1] * theta_1_distribution_array_array[:,1:2] * pragmatic_receiver_signal_not_1_h_prior
		
		print 'level %s pragmatic_receiver_signal_1_pre_normalized = \n%s' % (n, pragmatic_receiver_signal_1_pre_normalized)
		print 'level %s pragmatic_receiver_signal_not_1_pre_normalized = \n%s' % (n, pragmatic_receiver_signal_not_1_pre_normalized)
		
		pragmatic_receiver_signal_1_array = pragmatic_receiver_signal_1_pre_normalized / numpy.sum(pragmatic_receiver_signal_1_pre_normalized)
		pragmatic_receiver_signal_not_1_array = pragmatic_receiver_signal_not_1_pre_normalized / numpy.sum(pragmatic_receiver_signal_not_1_pre_normalized)
		
		pragmatic_receiver_memo[n] = numpy.asarray((pragmatic_receiver_signal_1_array, pragmatic_receiver_signal_not_1_array))

		pragmatic_sender_signal_1_h_array = numpy.sum(pragmatic_sender(n-1)[0] * theta_1_distribution_array_array[:,0:1], axis = 0)
		pragmatic_sender_signal_not_1_h_array = numpy.sum(pragmatic_sender(n-1)[1] * theta_1_distribution_array_array[:,1:2], axis = 0)

		pragmatic_sender_signal_1_fixed_theta_1_h_array = pragmatic_sender(n-1)[0][fixed_theta_1_num]
		pragmatic_sender_signal_not_1_fixed_theta_1_h_array = pragmatic_sender(n-1)[1][fixed_theta_1_num]

		pragmatic_receiver_signal_1_h_array = numpy.sum(pragmatic_receiver_memo[n][0], axis = 0)
		pragmatic_receiver_signal_1_h_array_densities = pragmatic_receiver_signal_1_h_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_1_h_array = \n%s' % pragmatic_receiver_signal_1_h_array
	
		pragmatic_receiver_signal_not_1_h_array = numpy.sum(pragmatic_receiver_memo[n][1], axis = 0)
		pragmatic_receiver_signal_not_1_h_array_densities = pragmatic_receiver_signal_not_1_h_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_not_1_h_array = \n%s' % pragmatic_receiver_signal_not_1_h_array
		
		pragmatic_receiver_signal_1_theta_1_array = numpy.sum(pragmatic_receiver_memo[n][0], axis = 1)
		pragmatic_receiver_signal_1_theta_1_array_densities = pragmatic_receiver_signal_1_theta_1_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_1_theta_1_array = \n%s' % pragmatic_receiver_signal_1_theta_1_array
	
		pragmatic_receiver_signal_not_1_theta_1_array = numpy.sum(pragmatic_receiver_memo[n][1], axis = 1)
		pragmatic_receiver_signal_not_1_theta_1_array_densities = pragmatic_receiver_signal_not_1_theta_1_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_not_1_theta_1_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_array
	
		fig, ax = pyplot.subplots(1, 2, figsize = (12,5))
	
		pyplot.subplot(1, 2, 1)
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_h_array, lw = 2, color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_h_array, lw = 2, color = 'c')
	
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'c')
	
		line = pyplot.plot(array_0, pragmatic_sender(n-1)[0][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender(n-1)[1][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'c')
	
		pyplot.subplot(1, 2, 2)
		line = pyplot.plot(array_0, pragmatic_receiver_signal_1_h_array_densities, lw = 2, color = 'b')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_1_theta_1_array_densities, lw = 2, linestyle = '--', color = 'b')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_not_1_h_array_densities, lw = 2, color = 'c')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_not_1_theta_1_array_densities, lw = 2, linestyle = '--', color = 'c')
	
		pyplot.subplot(1, 2, 1)
		pyplot.legend([r'$\sigma_{%s}(u_{1}|h)$' % (n-1), r'$\sigma_{%s}(\neg u_{1}|h)$' % (n-1), r'$\sigma_{%s}(u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(\neg u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2))], loc = 0, fontsize = 14)
	
		pyplot.subplot(1, 2, 2)
		pyplot.legend([r'$\rho_{%s}(h|u_{1})$' % n, r'$\rho_{%s}(\theta_{1}|u_{1})$' % n, r'$\rho_{%s}(h|\neg u_{1})$' % n, r'$\rho_{%s}(\theta_{1}|\neg u_{1})$' % n, ], loc = 0, fontsize = 14)

		fig.text(0, 0, r'$Lassiter\ and\ Goodman\ One\ Signal\ with\ Not\ Iterated\ Arbitrary$' + '\n', fontsize = 10)
	
		fig.text(0, 0, r'$\lambda = %s, C(u_{1}) \approx %s, C(\neg u_{1}) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{1}) \approx %s, C(\neg u_{1}) = %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + '\n' + r'$theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{1}) \approx %s, C(\neg u_{1}) = %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), num_states, theta_distribution_type) + '\n' + r'theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)
	
		# pyplot.savefig('Lassiter and Goodman Model One Signal with Not Iterated Arbitrary Normal Distribution.pdf')
		# pyplot.savefig('Lassiter and Goodman Model One Signal with Not Iterated Arbitrary Beta Distribution.pdf')
		# pyplot.savefig('Lassiter and Goodman Model One Signal with Not Iterated Arbitrary Uniform Distribution.pdf')
	
		pyplot.show()
		pyplot.close()

	return pragmatic_receiver_memo[n]

###################

# Here we have the settings for a level 0 receiver decoding probabilities, given a fixed
# theta. This forms the common basis for both Lassiter and Goodman's original model and
# our modified model.

cost = 2.
cost_of_not = 2./3.
choice_parameter =1.
lower_bound_list = [-4.]
upper_bound_list = [4.]
num_states = 80

mu = 0.
sigma = 1.
state_distribution = norm(mu,sigma)

# alpha_parameter = 1.
# beta_parameter = 9.
# location_parameter = lower_bound
# scale_parameter = upper_bound - lower_bound
# state_distribution = beta(alpha_parameter, beta_parameter, loc = location_parameter, scale = scale_parameter)

# state_distribution = uniform(lower_bound, upper_bound - lower_bound)

fixed_theta_1_num = numpy.int(numpy.ceil(num_states*(8./12.)))

theta_distribution_type = 'normal'
# theta_distribution_relation = 'True'
theta_posterior_source = 'signal_specific'

# pragmatic_sender_type = 'h_sensitive'
# pragmatic_sender_type = 'h_and_theta_sensitive'
# pragmatic_sender_type = 'h_sensitive_version_3'
pragmatic_sender_type = 'modified h sensitive'

if pragmatic_sender_type == 'modified h sensitive':
	weighting_sigma = 1.

if pragmatic_sender_type == 'modified h sensitive':
	weighting_array = numpy.empty((0, num_states))
	for h_num in numpy.arange(num_states):
		weighting_array = numpy.insert(weighting_array, h_num, truncnorm.pdf(array_0, lower_bound, upper_bound, loc = array_0[h_num], scale = weighting_sigma), axis = 0)
	print 'weighting_array = %s' % weighting_array
	
max_level = 10

#########################

array_0_array = numpy.empty((0, num_states))

bounds_pragmatic_receiver_h_array_densities = numpy.empty([0, max_level/2 + 1, 2, num_states])
for bound_num in numpy.arange(len(lower_bound_list)):
	lower_bound = lower_bound_list[bound_num]
	upper_bound = upper_bound_list[bound_num]

	if theta_distribution_type == 'normal':
		theta_1_distribution = norm(mu+2., sigma)
	elif theta_distribution_type == 'Beta':
		theta_1_distribution = beta(1, 9, loc = lower_bound, scale = upper_bound - lower_bound)
	elif theta_distribution_type == 'uniform':
		theta_1_distribution = uniform(lower_bound, upper_bound - lower_bound)

	array_0 = numpy.flipud(numpy.linspace(upper_bound, lower_bound, num_states, endpoint = False)) - ((numpy.flipud(numpy.linspace(upper_bound, lower_bound, num_states, endpoint = False)) - numpy.linspace(lower_bound, upper_bound, num_states, endpoint = False))/2)
	print 'array_0 = %s' % array_0
	
	array_0_array = numpy.insert(array_0_array, bound_num, array_0, axis = 0)
	
	theta_1_distribution_array = theta_1_distribution.pdf(array_0)
	theta_1_distribution_array = theta_1_distribution_array/numpy.sum(theta_1_distribution_array)
	initial_theta_1_distribution_array = theta_1_distribution_array
	initial_theta_1_distribution_array = numpy.reshape(initial_theta_1_distribution_array, [num_states, 1])

	print 'initial_theta_1_distribution_array = \n%s' % initial_theta_1_distribution_array
	print numpy.sum(initial_theta_1_distribution_array)

	plot_theta_1_distribution_array = initial_theta_1_distribution_array

	fig, ax = pyplot.subplots(1,1)
	pyplot.plot(array_0, plot_theta_1_distribution_array)
	pyplot.show()
	
	if pragmatic_sender_type == 'modified h sensitive':
		weighting_array = numpy.empty((0, num_states))
		for h_num in numpy.arange(num_states):
			weighting_array = numpy.insert(weighting_array, h_num, truncnorm.pdf(array_0, lower_bound, upper_bound, loc = array_0[h_num], scale = weighting_sigma), axis = 0)
	
	receiver_0_signal_1_array = numpy.empty([0, num_states])
	for theta_num in range(num_states):
		temp_signal_1_fixed_theta_array = numpy.empty(0)
		for h_num in range(num_states):
			value = receiver_0_signal_1(array_0[h_num], array_0[theta_num])
			temp_signal_1_fixed_theta_array = numpy.append(temp_signal_1_fixed_theta_array, value)
		receiver_0_signal_1_array = numpy.insert(receiver_0_signal_1_array, theta_num, temp_signal_1_fixed_theta_array, axis = 0)
	receiver_0_signal_1_array = receiver_0_signal_1_array / numpy.tile(numpy.sum(receiver_0_signal_1_array, axis = 1), (num_states, 1)).T
	receiver_0_signal_1_array = receiver_0_signal_1_array * initial_theta_1_distribution_array

	receiver_0_signal_not_1_array = numpy.empty([0, num_states])
	for theta_num in range(num_states):
		temp_signal_not_1_fixed_theta_array = numpy.empty(0)
		for h_num in range(num_states):
			value = receiver_0_signal_not_1(array_0[h_num], array_0[theta_num])
			temp_signal_not_1_fixed_theta_array = numpy.append(temp_signal_not_1_fixed_theta_array, value)
		receiver_0_signal_not_1_array = numpy.insert(receiver_0_signal_not_1_array, theta_num, temp_signal_not_1_fixed_theta_array, axis = 0)
	receiver_0_signal_not_1_array = receiver_0_signal_not_1_array / numpy.tile(numpy.sum(receiver_0_signal_not_1_array, axis = 1), (num_states, 1)).T
	receiver_0_signal_not_1_array = receiver_0_signal_not_1_array * initial_theta_1_distribution_array

	pragmatic_receiver_memo = {}
	pragmatic_sender_memo = {}
	pragmatic_receiver_memo[0] = numpy.asarray((receiver_0_signal_1_array, receiver_0_signal_not_1_array))

	pragmatic_receiver(max_level)
	levels_pragmatic_receiver_h_array_densities = numpy.empty([0, 2, num_states])
	for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):
		levels_pragmatic_receiver_h_array_densities = numpy.insert(levels_pragmatic_receiver_h_array_densities, pragmatic_receiver_level/2, numpy.sum(pragmatic_receiver_memo[pragmatic_receiver_level], axis = 1)/((upper_bound - lower_bound)/num_states), axis =0)
	bounds_pragmatic_receiver_h_array_densities = numpy.insert(bounds_pragmatic_receiver_h_array_densities, bound_num, levels_pragmatic_receiver_h_array_densities, axis = 0)

for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):

	fig, ax = pyplot.subplots(1, 1, figsize = (12,5))
	color_list = ['b', 'r', 'k']
	linestyle_list = ['-', '--', ':']
	pyplot.subplot(1, 1, 1)
	pyplot.grid(True)

	for bound_num in numpy.arange(len(lower_bound_list)):
		for signal in numpy.arange(2):
			pyplot.plot(array_0_array[bound_num], bounds_pragmatic_receiver_h_array_densities[bound_num, pragmatic_receiver_level/2, signal], color = color_list[signal], linestyle = linestyle_list[bound_num], lw = 3)
		pyplot.plot(array_0_array[bound_num], state_distribution.pdf(array_0_array[bound_num]), color = 'k', linestyle = '-', lw = 2)

	real_legend = pyplot.legend([r'$\rho_{%s}(h|u_{1})$' % pragmatic_receiver_level, r'$\rho_{%s}(h|\neg u_{1})$' % pragmatic_receiver_level, r'$\phi(h)$'], loc = 'lower left', bbox_to_anchor = (.75, .5))
	for legobj in real_legend.legendHandles:
		legobj.set_linewidth(1.5)

	ax = pyplot.gca().add_artist(real_legend)
	dummy_line_0, = pyplot.plot([], label = r'$a =-3, b =3$', color = 'gray', lw = 2, linestyle = '-')
	dummy_line_1, = pyplot.plot([], label = r'$a =-8, b =8$', color = 'gray', lw = 2, linestyle = '--')
	dummy_line_2, = pyplot.plot([], label = r'$a =-1, b =4$', color = 'gray', lw = 2, linestyle = ':')
	pyplot.legend(handles = [dummy_line_0, dummy_line_1, dummy_line_2], loc = 'lower right', bbox_to_anchor = (.25, .5), fontsize = 12)

	fig.text(0, 0, r'$Lassiter\ and\ Goodman\ One\ Signal\ with\ Not\ Iterated\ Arbitrary\ Variable\ Lower\ Bound\ Upper\ Bound$' + '\n', fontsize = 10)

	fig.text(0, 0, r'$\lambda = %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)
# 	fig.text(0, 0, r'$\lambda = %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)
# 	fig.text(0, 0, r'$\lambda = %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), num_states, theta_distribution_type) + r'theta\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, pragmatic_sender_type), fontsize = 10)


	pyplot.show()
	pyplot.close()
