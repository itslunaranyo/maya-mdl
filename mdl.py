# ==============================================================================

import maya.cmds as cmds
import qmdl.mdl as qmdl
import math, os, struct, imgbmp, shutil

projDir = "c:/projects/lunsp2/"
fileDestination = "c:/games/quake/lunsp2_dev/progs/"

def JustCopyADamnFile(src, dst):
	if not os.path.isfile(src):
		return
	srcdir = src.rsplit("/",1)[0] + "/"
	dstdir = dst.rsplit("/",1)[0] + "/"
	if not os.path.exists(dstdir):
		os.makedirs(dstdir)
	shutil.copy(src,dst)

# convert floats to strings and strip off the scientific notation from FPE
def niceroundstr( num ):
	out = str(round(num,5))
	if out.find("e") > 0:
		out = "0.0"
	return out

# maya is stupid and python is awesome and terrible
def vertexNormal( vert ):
	vfn = map(sum,zip(*zip(*[iter(cmds.polyNormalPerVertex( vert, q=1, xyz=1 ))]*3)))
	l = math.sqrt(vfn[0]*vfn[0] + vfn[1]*vfn[1] + vfn[2]*vfn[2])
	vfn = map(lambda x: x/l, vfn)
	return vfn

def crossProduct(a, b):
	return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]

def magnitude(vec):
	return math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])

# ================
	
class qcModel:
	def __init__(self):
		self.model = qmdl.Mdl()
		self.setDefaults()

	def __init__(self, scriptFile):
		self.model = qmdl.Mdl()
		self.setDefaults()
		self.exportScript(scriptFile)

	def setDefaults(self):
		self.name = 'model'
		self.meshes = ['mesh']
		self.base = None
		self.forwardnode = None
		self.skins = ['default']
		self.skindatas = []
		self.origin = '0 0 24'
		self.flags = 0
		self.skinpadding = 0
		self.scale = 1

		self.mdlProjDir = projDir
		self.mdlMayaBinDir = projDir+"models/"
		self.mdlSkinDir = projDir+"skins/"
		self.mdlFileDestination = fileDestination		
		
		self.shady = False
		self.totalFrames = 0
		self.uvbase = False
		
		self.tris = []
		self.verts = []
		self.frames = []
		self.framesmins = [999,999,999]
		self.framesmaxs = [-999,-999,-999]
		self.totalsize = 0
		
		self.framelist = []
	
	def crossTri(self, tri):
		a = self.verts[tri.vertices[0]]
		b = self.verts[tri.vertices[1]]
		c = self.verts[tri.vertices[2]]
		
		ab = [ 0, b.u-a.u, b.v-a.v ]
		cb = [ 0, b.u-c.u, b.v-c.v ]
		
		return crossProduct(ab, cb)	
	
	def SetFlags(self, fl):
		self.flags = int(fl)
		if (self.flags & 1):
			print("Setting Rocket Trail flag")
		if (self.flags & 2):
			print("Setting Smoke Trail flag")
		if (self.flags & 4):
			print("Setting Blood Trail flag")
		if (self.flags & 8):
			print("Setting Rotate flag")
		if (self.flags & 16):
			print("Setting Hellknight Trail flag")
		if (self.flags & 32):
			print("Setting Short Blood Trail flag")
		if (self.flags & 64):
			print("Setting Scrag Trail flag")
		if (self.flags & 128):
			print("Setting Vore Trail flag")
		if (self.flags & 0x4000):
			print("Setting Alphatest flag")
	
	# convert to quakespace: +x forward, +y left, +z up
	# I have been making all models in maya facing -z :|
	# TODO: make this switchable with a $command
	def quakeSpace(self, vector):
		return [-vector[2], -vector[0], vector[1]]
	
	# create a bunch of tuples pairing frame numbers to frame names to simplify formatFrames()
	def frameNames(self, frameStrings):
		mdlFrames = []
		frameTimes = []
		
		for seg in frameStrings[1:]:
			if len(seg) == 0:
				continue
			if seg.count("-") == 0:
				frameTimes.append(int(seg))
			else:
				start, end = seg.split("-")
				frameTimes.extend( range( int(start), int(end) + 1 ) )
		
		i=1
		for frame in frameTimes:
			mdlFrames.append([frame, frameStrings[0] + str(i)])
			i+=1
		
		return mdlFrames
	
	def frameForward(self, fwd):
		if (self.forwardnode is None):
			return 0, 0
		
		abs = cmds.xform(self.forwardnode, q=1, ws=1, t=1)[2] * self.scale * -1		# z movement
		rel = abs - fwd
		return abs, rel
	
	# take the vert positions on a model posed at current time and save them out
	# returns a 'frame' of the format [ ([vX, vY, vZ], [nX, nY, nZ]), ... ]
	def parseFrame( self, meshlist, s=None ):
		frame = []
		
		for m in meshlist:
			for vert in cmds.ls( cmds.polyListComponentConversion(m, tv=1), fl=1):
				vertOrg = cmds.xform(vert, q=1, ws=1, t=1)
				normal = vertexNormal(vert)
				if s >= 1:
					# smoke processing
					vertOrg = map(lambda (x, y): x + y * s, zip(vertOrg, normal))
					vertOrg[0] += math.cos( vertOrg[1] * 5 ) * s * 0.5
					vertOrg[1] += math.sin( vertOrg[1] * 5 ) * s * 0.5
					vertOrg[2] += s * s * 0.2
				
				#vertOrg = [-vertOrg[2], -vertOrg[0], vertOrg[1]]
				#normal = [-normal[2], -normal[0], normal[1]]
				vertOrg = self.quakeSpace(vertOrg)
				normal = self.quakeSpace(normal)
				
				# scale, then shift (so feet stay on the ground) to origin-relative coordinates
				vertOrg = map(lambda (x,y): x*self.scale-y, zip(vertOrg,self.origin))
				
				frame.append( (vertOrg, normal) )
			
				for k in xrange(3):
					self.framesmins[k] = min(self.framesmins[k],vertOrg[k])
					self.framesmaxs[k] = max(self.framesmaxs[k],vertOrg[k])
		return frame

	# parse the set of frames specified in frameNames
	# appends to self.frames the format (framelabel, frame, framefwd)
	def parseFrames( self, frameNames ):
		frames = 0
		fwd = 0			# total amount model has moved off the origin on mdl-x/mb-z, ROUNDED
		framefwd = 0	# difference in fwd from last frame, ROUNDED
		fshift = 0.0	# sub-unit amount the model has to be shifted if the forwards are scaled
		f = None
		
		for frame in frameNames:
			frameTime, frameLabel = frame
			oldorg = self.origin[:]
			
			# TODO: scaling for shade frames
			cmds.currentTime( frameTime )
			if self.shady:
				for s in range(1,6):
					f = self.parseFrame( self.meshes, s )
					self.frames.append( (frameLabel + "_" + str(s), f ) )
				frames += 5
			else:
				fwd, framefwd = self.frameForward(fwd)
			#	fshift = framefwd - round(framefwd)
		#		self.origin[0] = fshift + oldorg[0]		# bad hack
				f = self.parseFrame( self.meshes )
				# FIXME: framefwd sometimes shows a cumulative 1 unit error on the last frame of a loop
		#		self.frames.append( ( frameLabel, f, round(framefwd) ) )
				self.frames.append( ( frameLabel, f, round(framefwd, 2) ) )
				frames += 1
		self.origin = oldorg
		return frames

	# parse the set of frames specified in frameNames as a frameGroup
	# appends to self.frames the format (framelabel, [frames])
	def parseFrameGroup( self, groupName, frameNames ):
		fg = []
	
		for frame in frameNames:
			frameTime, frameLabel = frame
			cmds.currentTime( frameTime )
			
			if self.shady:
				for s in range(1,6):
					fg.append( (frameLabel, self.parseFrame( self.meshes, s ), 0 ) )
			else:
				fg.append( ( frameLabel, self.parseFrame( self.meshes ), 0 ) )
		
		self.frames.append( (groupName, fg) )
		return 1

	def storeFaceVertIndices(self):
		print "cacheing face->vert indices ..."
		offset = 0
		for m in self.meshes:
			for poly in cmds.ls( cmds.polyListComponentConversion(m, tf=1), fl=1):
				polyVtxFaces = cmds.ls( cmds.polyListComponentConversion(poly, ff=1, tvf=1), fl=1)
				if len(polyVtxFaces) is not 3:
					raise Exception("You didn't triangulate your model!")
					
				if self.uvbase:
					triVerts = []
					for vtf in polyVtxFaces:
						uv = cmds.ls( cmds.polyListComponentConversion(vtf, fvf=1, tuv=1), fl=1)[0]
						triVerts.append( offset + int( uv.partition("[")[2][:-1] ) )
				else:
					triVerts = map( lambda x: offset + int( x.partition("[")[2].partition("]")[0] ), polyVtxFaces )
				
				triVerts.reverse()
				newTri = qmdl.Mdl.Triangle()
				newTri.vertices = tuple( triVerts )
				self.tris.append( newTri )
			offset += cmds.polyEvaluate(m,uv=1)
		
	def parseBaseUV(self):
		print "parsing base frame from UVs ..."
		self.verts = []
		self.UVsToVerts = []
		offset = 0
		
		for m in self.meshes:
			for uv in cmds.ls( cmds.polyListComponentConversion(m, tuv=1), fl=1):
				self.UVsToVerts.append(offset + int(cmds.polyListComponentConversion(uv, fuv=1, tv=1)[0].partition("[")[2][:-1]))
				uvOrg = cmds.polyEditUV(uv, q=1)
				
				newVert = qmdl.Mdl.Vertex()
				newVert.u, newVert.v = uvOrg
				self.verts.append( newVert )
				newVert.u = math.floor(newVert.u * self.model.skinwidth + 0.5)
				newVert.v = math.floor((1-newVert.v) * self.model.skinheight + 0.5)
			offset += cmds.polyEvaluate(m,v=1)

		self.storeFaceVertIndices()
		
		for tri in self.tris:
			tri.backface = qmdl.FACE_FRONT
			cx = self.crossTri(tri)
			self.totalsize += magnitude(cx) * 0.5
			
		self.model.size = self.totalsize * self.model.skinwidth * self.model.skinheight
			
	def parseBase(self):
		print "parsing base frame ..."
		cmds.currentTime( 1 )
		self.verts = []
		basex = 0
		basemins = [999,999]
		basemaxs = [-999,-999]
		
		# convert base pose to s & t
		for vert in cmds.ls( cmds.polyListComponentConversion(self.base, tv=1), fl=1):
			vertOrg = cmds.xform(vert, q=1, ws=1, t=1)
			# convert to quakespace: +x forward, +y left, +z up
			vertOrg = [-vertOrg[2], -vertOrg[0], -vertOrg[1]]
			
			newVert = qmdl.Mdl.Vertex()
			newVert.u, newVert.v = vertOrg[1], vertOrg[2]
			
			basemins[0] = min(basemins[0],newVert.u)
			basemaxs[0] = max(basemaxs[0],newVert.u)
			basemins[1] = min(basemins[1],newVert.v)
			basemaxs[1] = max(basemaxs[1],newVert.v)

			self.verts.append( newVert )

		self.storeFaceVertIndices()
		
		# calculate facings & seams
		for tri in self.tris:
			cx = self.crossTri(tri)
			
			if cx[0] > 0:
				tri.backface = qmdl.FACE_BACK
				self.verts[tri.vertices[0]].onseam = -1
				self.verts[tri.vertices[1]].onseam = -1
				self.verts[tri.vertices[2]].onseam = -1
			else:
				tri.backface = qmdl.FACE_FRONT
			
			# calculate size here too even though only software quake uses it
			self.totalsize += magnitude(cx) * 0.5
		
		for tri in self.tris:
			a = self.verts[tri.vertices[0]]
			b = self.verts[tri.vertices[1]]
			c = self.verts[tri.vertices[2]]
			if tri.backface == qmdl.FACE_FRONT:
				a.onseam = abs(a.onseam)
				b.onseam = abs(b.onseam)
				c.onseam = abs(c.onseam)
		
		basescale = [basemaxs[0]-basemins[0], basemaxs[1]-basemins[1]]
		for vert in self.verts:
			vert.u = (vert.u - basemins[0])
			vert.v -= basemins[1]
		
		for vert in self.verts:
			vert.u /= basescale[0]
			vert.v /= basescale[1]
		
			vert.u = self.skinpadding + math.floor(vert.u * (self.model.skinwidth * 0.5 - 2 * self.skinpadding) + 0.5)
			vert.v = self.skinpadding + math.floor(vert.v * (self.model.skinheight - 2 * self.skinpadding) + 0.5)
			
			if vert.onseam == -1:
				vert.onseam = 0
				vert.u += self.model.skinwidth / 2
				
		self.model.size = self.totalsize * self.model.skinwidth * self.model.skinheight / (basescale[0] * basescale[1])
	
	def calculateCoordsFrame( self, frame ):
		name, vs, fwd = frame
		newFrame = qmdl.Mdl.Frame()
		
		if self.uvbase:
			uvvs = map(lambda x: vs[self.UVsToVerts[x]], xrange(len(self.verts)))
			vs = uvvs

		for v in vs:
			c = qmdl.Mdl.Coord()
			point, normal = v
			newPos = [0,0,0]
			for i in xrange(3):
				newPos[i] = math.floor( (point[i] - self.model.origin[i]) / self.model.scale[i] )
			c.position = tuple(newPos)
			c.encode(normal)
			newFrame.vertices.append(c)
		newFrame.name = name
		newFrame.calculate_bounds()
		
		return newFrame
	
	def calculateCoordsFrameGroup( self, frameGroup ):
		name, frames = frameGroup
		newFrameGroup = qmdl.Mdl.FrameGroup()
		delay = 0.1
		for frame in frames:
			newFrameGroup.duration.append(delay)
			delay += 0.1
			newFrameGroup.frames.append(self.calculateCoordsFrame(frame))
		
		return newFrameGroup		
		
	def calculateCoords( self ):
		print "compressing vertex coordinates ..."
		frames = []
		for frame in self.frames:
			# TODO: change qmdl, add raw coordinates?
			if len(frame) == 2: # need better way to tell the difference between frame and framegroup but type checking is Bad
				frames.append( self.calculateCoordsFrameGroup(frame) )
			else:
				frames.append( self.calculateCoordsFrame(frame) )
			
		return frames
	
	def writeMDL(self):
		boundleg = [0,0,0]
		scl = [0,0,0]

		for i in xrange(3):
			scl[i] = (self.framesmaxs[i] - self.framesmins[i]) / 255.9;
			boundleg[i] = max(abs(self.framesmins[i]), abs(self.framesmaxs[i]) )
		self.model.scale = tuple(scl)
		self.model.origin = tuple(self.framesmins)

		self.model.boundingradius = magnitude(boundleg)

		self.model.vertices = self.verts
		self.model.triangles = self.tris
		self.model.frames = self.calculateCoords()
		print("Total verts: " + str(len(self.verts)))
		print("Total tris: " + str(len(self.tris)))
		print("Total frames: " + str(len(self.model.frames)))

		self.model.flags = self.flags
		self.model.synctype = 1

		fdest = self.mdlFileDestination + self.name + ".mdl"
		JustCopyADamnFile(fdest, fdest + ".bak")
		mdlDest = open( fdest, "w+" )
		
		print "writing " + fdest
		self.model.write(mdlDest)
		os.utime(fdest, None)
		mdlDest.close()
		print (str(self.totalFrames) + " frames written to " + self.mdlFileDestination + self.name + ".mdl")

	def readSkins(self):
		print "reading skins ..."
		self.model.skinwidth, self.model.skinheight = 0, 0
		self.model.skins = []

		# skins have to be windows 8-bit non-compressed non-row-flipped .BMP files, because
		# this code is astoundingly lazy, because I can't be bothered to mess with image
		# libraries. copies raw binary data right out of the bitmap with no conversion.
		# if you save the file properly, this works 100% of the time, and if you don't you
		# get utter garbage with no warnings glhf
		for skin in self.skins:
			skfile = open(self.mdlSkinDir + skin + ".bmp","r")
			skimg = imgbmp.imgBMP()
			skimg.load(skfile)
			if self.model.skinwidth == 0 or (self.model.skinwidth == skimg.width and self.model.skinheight == skimg.height):
				newSkin = qmdl.Mdl.Skin(skimg.width,skimg.height)
				self.model.skinwidth = skimg.width
				self.model.skinheight = skimg.height
				#print skimg.pdata
				newSkin.pixels = skimg.pdata
				self.model.skins.append(newSkin)
	
	def printQCString(self):
		# framestr is the block of $frame defs that go at the top of the .qc file
		framestr = ""
		# funcstr is the block of [frame,nextthink] frame macros
		funcstr = ""
		
		f = 0
		for framegroupsize in self.framelist:
			framestr += "$frame "
			for i in range(framegroupsize):
				framename = self.frames[f][0]
				framefwd = self.frames[f][2]
				framestr += framename + " "
				
				fwdstr = ""
				if (self.forwardnode is not None):
					if (framefwd < 0):
						framefwd = abs(framefwd)
						fwdstr = "ai_forward(" + str(framefwd) + ");"
					elif (framefwd > 0):
						fwdstr = "ai_back(" + str(framefwd) + ");"
				
				funcstr += "void()	" + self.name + "_" + framename + " =	[ $" + framename + ",	" + self.name + "_"
				
				if (i+1 == framegroupsize):
					funcstr += self.frames[f+1-framegroupsize][0]
				else:
					funcstr += self.frames[f+1][0]
				
				funcstr += " ]	{ " + fwdstr + " }\n"
				f += 1
			
			framestr += "\n"
			funcstr += "\n"
		
		print("\n// ================================\n\n")
		print framestr
		print funcstr
	
	def fixPath(self, pth):
		p = pth[:]
		p = p.strip("\'\"").replace("\\","/")
		if p[-1] != "/":
			p += "/"
		return p
	
	# load a very qc-like .txt file of model export script commands
	# open the maya files it lists one by one
	# convert their posed models to frames in a .mdl
	def exportScript(self, scriptfilename):
		path1 = self.mdlMayaBinDir + scriptfilename
		print("opening", path1)
		scriptfile = open( path1 )
		scriptlines = list(scriptfile)
		scriptfile.close()
		
		pastHeader = False
		baseFrameWritten = False
		print ("parsing")
		# parse file for commands beginning with $
		linenum = 0
		for line in scriptlines:
			linenum += 1
			if line[0] != "$":
				continue
			
			tokens = line.strip("$\r\n").split(" ")
			cmd = tokens[0]
			ph_cmds = ["file", "anim", "animgroup", "set"]
			
			if (pastHeader == False):
				if cmd == "name":
					self.name = tokens[1]
				elif cmd == "mesh":
					self.meshes = tokens[1:]
				elif cmd == "basemesh":
					self.base = tokens[1]
				elif cmd == "forwardnode":
					self.forwardnode = tokens[1]
				elif cmd == "useuvs" or cmd == "useUVs":
					self.uvbase = True
				elif cmd == "scale":
					self.scale = float(tokens[1])
				elif cmd == "skins":
					self.skins = tokens[1:]
				elif cmd == "skinpadding":
					self.skinpadding = int(tokens[1])
				elif cmd == "shady":
					self.shady = True
				elif cmd == "origin":
					self.origin = map(float,tokens[1:4])
				elif cmd == "flags":
					self.SetFlags(tokens[1])
				elif cmd.lower() == "projdir":
					self.mdlProjDir = self.fixPath(tokens[1])
				elif cmd.lower() == "skindir":
					self.mdlSkinDir = self.fixPath(tokens[1])
				elif cmd.lower() == "filedestination":
					self.mdlFileDestination = self.fixPath(tokens[1])
				elif cmd in ph_cmds:
					pastHeader = True
				else:
					print("Bad command token '" + cmd + "' in header on line " + str(linenum))
					return
			
			# TODO: eliminate this requirement - run to EOF twice, do header commands in first pass and these second
			if (pastHeader == True):
				print("parsing anims")
				if cmd not in ph_cmds:
					print("Bad command token '" + cmd + "' outside header on line " + str(linenum))
					return
			
				if cmd == "file":
					currentFile = cmds.file(q=1,sn=1).lower()
					neededFile = (self.mdlProjDir + tokens[1]).lower()
					print(currentFile, "\n", neededFile)
					if currentFile != neededFile:
						print ("opening " + neededFile)
						try:
							mayaScene = cmds.file(self.mdlProjDir + tokens[1], f=1, options="v=0", o=1)
						except RuntimeError:
							raise Exception(neededFile + " could not be opened for some bullshit reason")
						
					if not baseFrameWritten:
						self.readSkins()		
						if self.uvbase:
							if self.base:
								raise Exception("export script can't have both basemesh and baseUVs")
							self.parseBaseUV()
						else:
							if not self.base:
								raise Exception("export script missing both basemesh and baseUVs")
							self.parseBase()
							
						baseFrameWritten = True
				
				elif cmd == "animgroup":
					added = self.parseFrameGroup( tokens[1], self.frameNames(tokens[1:]) )
					self.totalFrames += added
					self.framelist.append( added )
				elif cmd == "anim":
					added = self.parseFrames( self.frameNames(tokens[1:]) )
					self.totalFrames += added
					self.framelist.append( added )
				elif cmd == "set":
					print("setting " + tokens[1] + " to " + tokens[2])
					cmds.setAttr(tokens[1],float(tokens[2]))
				# framelist is a list of frame group sizes so printQCString knows where to loop without
				# guessing at frame name/number suffixes
				
		print (str(self.totalFrames) + " frames grabbed from scenes")
		if self.totalFrames > 255:
			raise Exception("that's too many frames! jesus!")
		
		self.writeMDL()
		
# ================
# modifications to QMDL:
# todo: this SHOULD be a new class that inherits from qmdl.Mdl

class CoordRaw:
	"""
	This class stores an uncompressed coordinate before being processed into
	a compressed mdl.
	
	"""
	def __init__(self):
		"""Initialise all the class members with valid values."""
		self.position, self.normal = (0, 0, 0), (0, 0, 0)

	def compress(self, origin, scale):
		c = Coord()
		
		cPos = [0,0,0]
		for i in xrange(3):
			cPos[i] = math.floor( (self.position[i] - origin[i]) / scale[i] )
		c.position = tuple(cPos)
		c.encode(self.normal)
		return c

qmdl.Mdl.CoordRaw = CoordRaw

def frame__init(self):
	"""
	Initialise all the members

	"""
	self.bbox_min = qmdl.Mdl.Coord()
	self.bbox_max = qmdl.Mdl.Coord()
	self.name = b"nameless"
	self.vertices = []
	self.verticesRaw = []	# lunaran change

def frame__calculate_bounds(self):
	"""
	Find the compressed bounds of this frame, if none was read from a file.
	"""
	
	boxmin = [999,999,999]
	boxmax = [-999,-999,-999]
	for c in self.vertices:
		for k in xrange(3):
			boxmin[k] = min(boxmin[k],c.position[k])
			boxmax[k] = max(boxmax[k],c.position[k])
	self.bbox_min.position = tuple(boxmin)
	self.bbox_max.position = tuple(boxmax)

def frame__calculate_bounds_raw(self):
	"""
	Return the (uncompressed) bounds of this frame's raw coordinates.
	"""
	mins = [999,999,999]
	maxs = [-999,-999,-999]
	for v in self.verticesRaw:
		for k in xrange(3):
			mins[k] = min(mins[k],v.position[k])
			maxs[k] = max(maxs[k],v.position[k])
	return mins, maxs

def frame__compress(self, origin, scale):
	"""
	Converts all the raw floating point coordinates in a frame to integerized
	Coord objects. Frames don't know about other frames, and only the mdl
	object itself knows the maximum extents of all frames, so the origin and
	scale parameters which determine how each frame is integerized come from
	outside as parameters.
	
	"""
	self.vertices = map(lambda x: x.compress(origin,scale), self.verticesRaw)
	self.calculate_bounds()

qmdl.Mdl.Frame.__init__ = frame__init
qmdl.Mdl.Frame.calculate_bounds = frame__calculate_bounds
qmdl.Mdl.Frame.calculate_bounds_raw = frame__calculate_bounds_raw
qmdl.Mdl.Frame.compress = frame__compress

def framegroup__compress(self, origin, scale):
	"""
	Compress all frames in the group.  Since Frames and FrameGroups are stored
	in the same flat list, FrameGroups must provide a compress() method with
	an identical signature.
	
	"""
	for frame in self.frames:
		frame.compress(origin, scale)
	self.calculate_bounds()

qmdl.Mdl.FrameGroup.compress = framegroup__compress

def calculate_bounds_raw(self):
	mins, maxs = [999,999,999],[-999,-999,-999]
	
	for frame in self.frames:
		fmin, fmax = frame.calculate_bounds_raw()
		for i in xrange(3):
			mins[i] = min(mins[i], fmin[i])
			maxs[i] = max(maxs[i], fmax[i])
	
	return mins, maxs			

def compress(self):
	boundleg = [0,0,0]
	scl = [0,0,0]

	mins, maxs = self.calculate_bounds_raw()
	
	for i in xrange(3):
		scl[i] = (maxs[i] - mins[i]) / 255.9;
		boundleg[i] = max(abs(mins[i]), abs(maxs[i]) )

	self.origin = tuple(self.framesmins)
	self.scale = tuple(scl)
	self.boundingradius = magnitude(boundleg)

	for frame in self.frames:
		frame.compress(self.origin, self.scale)

qmdl.Mdl.calculate_bounds_raw = calculate_bounds_raw
qmdl.Mdl.compress = compress

# ==============================================================================
# MESH IMPORT

# importMDL("d:/games/quake/lunsp2/oldprogs/knight.mdl")

import struct

class vertex:
	def __init__(self, (v1, v2, v3, normal)):
		self.v = [v1, v2, v3] # [v1, -v3, v2] # convert to maya's y-up
		self.normal = normal
	
	def resize(self,scale,translate):
		for i in range(3):
			self.v[i] = self.v[i] * scale[i] + translate[i]
	
	def objv(self):
		return "v " + str(self.v[0]) + " " + str(self.v[2]) + " " + str(-self.v[1]) + "\n"

class frame:
	def __init__(self, type, vmin, vmax, name):
		self.mins = vmin
		self.maxs = vmax
		self.name = name[0].split('\x00')[0]
		self.verts = []
	
	def vert(self, v):
		self.verts.append(v)
	
	def resize(self,scale,translate):
		for v in self.verts:
			v.resize(scale,translate)
	

def readFrom(file, pattern):
	return struct.unpack(pattern, file.read(struct.calcsize(pattern)))

def readVertexFrom(file):
	return vertex(struct.unpack("4B", file.read(4)))

def importMDL(infile, maxFrames=0):
	FH = open(infile, 'rb')

	# read header info
	ident,version = readFrom(FH, "ii")
	modelscale = readFrom(FH, "3f")
	modeltranslate = readFrom(FH, "3f")
	boundingradius = readFrom(FH, "f")
	eyeposition = readFrom(FH, "3f")
	skinvals = readFrom(FH, "3i")
	numVerts,numTris,numFrames = readFrom(FH, "iii")
	if maxFrames:
		if maxFrames < numFrames:
			numFrames = maxFrames
	synctype, flags, size = readFrom(FH, "iif")

	# skip the textures
	# TODO: we're assuming no skingroups
	skinbytes = skinvals[0] * (4 + skinvals[1]*skinvals[2])
	FH.seek(skinbytes, 1)

	# skip the STs
	stbytes = numVerts*12
	FH.seek(stbytes, 1)

	tris = []
	for i in range(numTris):
		tris.append(readFrom(FH, "4i"))

	frames = []
	for i in range(numFrames):
		fr = frame(readFrom(FH, "i"),readVertexFrom(FH),readVertexFrom(FH),readFrom(FH, "16s"))
		for j in range(numVerts):
			fr.vert(readVertexFrom(FH))

		fr.resize(modelscale,modeltranslate)
		frames.append(fr)

	FH.close()

	tempfile = infile + ".obj"

	k = 1
	for f in frames:
		OBJ = open(tempfile, 'w')
		for v in f.verts:
			OBJ.write(v.objv())
		
		for t in tris:
			OBJ.write("f " + str(t[1]+k) + " " + str(t[2]+k) + " " + str(t[3]+k) + " " + "\n")
		
		#k += numVerts

		OBJ.close()
		cmds.file(tempfile, i=1, type="OBJ", ra=1, rpr=f.name, options="mo=1", pr=1)
	
	os.remove(tempfile)




def uvs2flatverts():
	for obj in cmds.ls(sl=1,fl=1):
		uv = cmds.polyListComponentConversion(obj,fv=1,tuv=1)
		st = cmds.polyEditUV(uv,q=1)
		xyz = cmds.xform(obj,t=1,q=1)
		xyzuv = (st[0]*128, st[1]*128, xyz[2])
		cmds.xform(obj,r=0,ws=1,t=xyzuv)