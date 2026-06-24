import numpy

# Read data
data=numpy.loadtxt('output_000000.dat', usecols=(1,9,10))

# Find zero crossings
zc = numpy.where(numpy.diff(numpy.sign(data[:,2])))[0]

# Get time and coefficients
cd=-2*data[:,1]
cl=2*data[:,2]

# Get force and time before and after crossing
t0=data[zc,0]
t1=data[zc+1,0]

f0=data[zc,2]
f1=data[zc+1,2]

# Interpolate to find actual time of zero crossing
dt = -f0/(f1-f0)
t = t0 + dt*(t1 - t0)
T = t[2:] - t[:-2]

# Report average and fluctuating components of lift and drag
# Note that the period of the drag is twice that of the lift
for i in range(0,int(numpy.log2(len(zc)/2))):
   s = 2**i
   print("------------------------------------------")
   print("Average over", s, "periods")
   print("------------------------------------------")
   ss =  zc[-1] - zc[-2*s-1]

   fft = numpy.fft.fft(cd[zc[-2*s-1]:zc[-1]])
   print("Cd  =", round(  numpy.abs(fft[0  ]/ss),5),
         "Cd' =", round(2*numpy.abs(fft[2*s]/ss),5),
         "p2p =", round(4*numpy.abs(fft[2*s]/ss),5))

   fft = numpy.fft.fft(cl[zc[-2*s-1]:zc[-1]])
   print("Cl  =",round(  numpy.abs(fft[0]/ss),6),
         "Cl' =",round(2*numpy.abs(fft[s]/ss),5),
         "p2p =",round(4*numpy.abs(fft[s]/ss),5))

   print("Check = ",
         round(numpy.average(cd[zc[-2*s-1]:zc[-1]]),5),
         round(numpy.average(cl[zc[-2*s-1]:zc[-1]]),5))
   print("Period   = ",round(numpy.average(T[-s:]),5))
   print("Strouhal = ",round(1.0/numpy.average(T[-s:]),5))


print("==========================================")
print("Engelman & Jamnia 1990")
print("==========================================")
print("Cd  =", round(1.411,5),
      "Cd' =", round(0.0203/2,5),
      "p2p =", round(0.0203,5))
print("Cl  =",round(0.0,6),
      "  Cl' =",round(0.7267/2,5),
      "p2p =",round(0.7267,5))
print("Period   = ",round(5.80,5))
print("Strouhal = ",round(0.173,5))
print("==========================================")



