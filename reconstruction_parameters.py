class reconstruction_parameters:
	def __init__(self, smoothing, OmegaM0, bias, f):
		self.smoothing = smoothing
		self.OmegaM0 = OmegaM0
		self.bias = bias
		self.f = f
	def OmegaM(self, z):
		Ez2 = self.OmegaM0*(1+z)**3 + (1-self.OmegaM0)
		return self.OmegaM0 * (1+z)**3/Ez2
