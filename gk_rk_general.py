def cross_cor_coeff(delta_k_1,pk_1,delta_k_2,pk_2,L,N,pkxibins_spec):
	import numpy as np
	import correlation_general as cor
	

	k_avg,correlation_avg=cor.correlation(delta_k_2,delta_k_1,L,N,pkxibins_spec) # numerator

	print("correlation_avg shape=", np.shape(correlation_avg))
	print("pk_2 shape=",np.shape(pk_2))
	print("pk_1 shape=", np.shape(pk_1))

	corr_coeff_dis_rk=correlation_avg/(pk_2*pk_1)**(1./2.)

	corr_coeff_dis_gk=correlation_avg/pk_1

	return k_avg,corr_coeff_dis_gk,corr_coeff_dis_rk,correlation_avg


