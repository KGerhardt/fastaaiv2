import gzip

class agnostic_reader:
	def __init__(self, file):
		self.file = file
		
		self.dat = None
		self.is_gz = None
	
	def read(self):
		with open(self.file, 'rb') as fh:
			self.dat = fh.read()
			
		if len(self.dat) > 2:
			self.is_gz = self.dat[0:2] == b'\x1f\x8b'
		else:
			self.is_gz = None #This means the file is empty, or at least has less than 2 bytes of info.
			self.dat = None
			
		if self.is_gz is not None:
			if self.is_gz:
				self.dat = gzip.decompress(self.dat)
			
			try:
				self.dat = self.dat.decode(encoding = 'ascii')
			except:
				self.dat = self.dat.decode(encoding = 'utf-8')
		
		return self.dat
	
		