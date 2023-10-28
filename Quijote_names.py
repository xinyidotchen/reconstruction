class Quijote_names:
	def __init__(self,category,subcategory,sim_num,grid_size):
		self.category=category
		self.subcategory=subcategory
		self.sim_num=sim_num
		self.grid_size=grid_size
	def whoami(self):
		print("This is Quijote %s %s %i"%(self.category,self.subcategory,self.sim_num,self.grid_size,))
	def getFileName(self):
		return "Quijote_%s_%s_%06i"%(self.category,self.subcategory,self.sim_num,)
	def get_bf_rec_DirName(self,redshift):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/bf_rec"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,)
	def get_bf_rec_Halos_DirName(self,redshift,mass_lim_low,mass_lim_up):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/bf_rec/%.2f-%.2f"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,mass_lim_low,mass_lim_up)

	def get_ini_DirName(self):
		return "Quijote/Snapshots/%s/%i/grid%i/ini"%(self.subcategory,self.sim_num,self.grid_size,)
	def get_hada_DirName(self,redshift,smooth_scale):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/hada/sm%i"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,smooth_scale)
	def get_hada_Halos_DirName(self,redshift,mass_lim_low,mass_lim_up,smooth_scale):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/hada/%.2f-%.2f/sm%i"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,mass_lim_low,mass_lim_up,smooth_scale)

	def getSimPath(self):
		return "/gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i"%(self.category,self.subcategory,self.sim_num,)
	def get_hada_header(self,redshift,space,smooth_scale,Cani):
		return '''
			# Simulation path : /gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i
			# Simulation redshift : %.1f
			# Grid size : %i
			# Space: %s
			# Hada method
			# Eff smoothing: %.1f Mpc/h
			# Cani: %.1f\n
			'''%(self.category,self.subcategory,self.sim_num,redshift,self.grid_size,space,smooth_scale,Cani,)
	def get_hada_Halos_header(self,redshift,space,mass_lim_low,mass_lim_up,smooth_scale,Cani):
		return '''
			# Simulation path : /gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i
			# Simulation redshift : %.1f
			# Grid size : %i
			# Space: %s
			# mass bin: %.2f-%.2f
			# Hada method
			# Eff smoothing: %.1f Mpc/h
			# Cani: %.1f\n
			'''%(self.category,self.subcategory,self.sim_num,redshift,self.grid_size,space,mass_lim_low,mass_lim_up,smooth_scale,Cani,)


	def get_standard_DirName(self,redshift,smooth_scale):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/standard/sm%i"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,smooth_scale)
	def get_standard_Halos_DirName(self,redshift,mass_lim_low,mass_lim_up,smooth_scale):
		return "Quijote/%s/%s/%i/grid%i/z%.1f/standard/%.2f-%.2f/sm%i"%(self.category,self.subcategory,self.sim_num,self.grid_size,redshift,mass_lim_low,mass_lim_up,smooth_scale)

	def get_standard_header(self,redshift,space,smooth_scale):
		return '''
			# Simulation path : /gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i
			# Simulation redshift : %.1f
			# Grid size : %i
			# Space: %s
			# Standard method
			# Smoothing: %.1f Mpc/h\n
			'''%(self.category,self.subcategory,self.sim_num,redshift,self.grid_size,space,smooth_scale,)
	def get_standard_Halos_header(self,redshift,space,mass_lim_low,mass_lim_up,smooth_scale):
		return '''
			# Simulation path : /gpfs/loomis/scratch60/padmanabhan/np274/simulations/quijote/%s/%s/%i
			# Simulation redshift : %.1f
			# Grid size : %i
			# Space: %s
			# mass bin: %.2f-%.2f
			# Standard method
			# Smoothing: %.1f Mpc/h\n
			'''%(self.category,self.subcategory,self.sim_num,redshift,self.grid_size,space,mass_lim_low,mass_lim_up,smooth_scale,)


