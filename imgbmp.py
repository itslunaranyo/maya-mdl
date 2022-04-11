import struct

# quick and dirty palettized bitmap ripper

def _struct_read(s, file):
    return s.unpack(file.read(s.size))

class imgBMP:
	def __init__(self):
		self.width = 0
		self.height = 0
		self.pdata = ''
		
	def load(self, bmpfile):
		fileheaderstruct = struct.Struct("<HiHHi")
		infoheaderstruct = struct.Struct("<illhhiillii")
		palettestruct = struct.Struct("<256l")
		fheader = _struct_read(fileheaderstruct, bmpfile)
		iheader = _struct_read(infoheaderstruct, bmpfile)
		palette = _struct_read(palettestruct, bmpfile)

		self.width = iheader[1]
		self.height = iheader[2]
		self.pdata = b''
		
		bmpfile.seek(self.width * (self.height - 1), 1)
		
		for i in range(self.height):
			self.pdata += bmpfile.read(self.width)
			bmpfile.seek(-self.width*2, 1)
	
	def crop(self, top, left, bottom, right):
		subbmp = imgBMP()
		if top < 0 or bottom > self.height or left < 0 or right > self.width:
			print("bad bounds on crop")
			return None
		
		subbmp.width = right - left
		subbmp.height = bottom - top
		subbmp.pdata = b''
		
		for row in range(top, top+subbmp.height):
			subbmp.pdata += self.pdata[row * self.width + left : row*self.width + right]
		
		return subbmp
	
	def writeLMP(self, file):
		lmpstruct = struct.Struct("<ll"+str(self.width*self.height)+"s")
		lmpdat = lmpstruct.pack(self.width, self.height, self.pdata) 
		with open(file, "wb+") as lmpfile:
			lmpfile.write(lmpdat)