import matplotlib.pyplot as plt

xRange = range( 0, 8192 )

yBeforePhaserot = list()
yAfterPhaserot = list()

def ReadFileFloats( filename: str ) -> list:
	global xRange

	with open( filename ) as txtfile:
		floats = [ float( curline ) for curline in txtfile ]
		floats = floats[ : len( xRange ) ]

	return floats

y_rot = ReadFileFloats( "graph_before_phaserotation.txt" )
y_lpf = ReadFileFloats( "graph_before_lowpassfilter.txt" )
y_fd1 = ReadFileFloats( "graph_before_freqdoubling1.txt" )
y_end = ReadFileFloats( "graph_before_end.txt" )

#y_justdecoded = ReadFileFloats( "just_decoded.txt" )
#plt.plot( xRange, y_justdecoded, linewidth = 0.5, c = 'black', label='Just decoded' )

plt.plot( xRange, y_rot, linewidth = 0.5, c = 'darkred', label='1. Initial' )
plt.plot( xRange, y_lpf, linewidth = 0.5, c = 'orangered', label='2. After phaserotation' )
plt.plot( xRange, y_fd1, linewidth = 0.5, c = 'yellow', label='3. After low pass filter' )
plt.plot( xRange, y_end, linewidth = 0.5, c = 'forestgreen', label='4. After freq doubling' )
# '''

plt.margins( 0.0 )
plt.tight_layout()
plt.legend()
plt.show()
