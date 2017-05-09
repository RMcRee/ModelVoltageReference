# ModelVoltageReference
Model pink noise and temperature dependence of voltage references

Short discussion of my purpose in creating this here:
http://www.eevblog.com/forum/metrology/software-simulation-of-a-voltage-reference-noise-tempco-long-term-drift/

The idea is to simulate some of the more salient properties of voltage references that affect their performance,
namely noise and temperature.

For example, let's say that you build a circuit which can measure the difference between two voltage references and 
amplify that difference (e.g. an instrumentation amplifier). The simulation ought to give us some insight into the
noise and temperature dependence of the INAMP output. Furthermore, if we have a set of, say, four references, and know their starting values can we track their changes over temperature and long-term drift using just differential measurements?

Simulation shows that differential-only measurements *are* a viable way to accurately measure a set of references. After developing this, I googled "NIST least squares metrology" and came up with this 1974 paper which discusses measuring differences in gage blocks (coming up with an almost identical least squares equation):
http://emtoolbox.nist.gov/Publications/NBSIR74-587.asp
See especially appendix at end. Note that author talks about adding constraint to LSQ: (A+B)=K on page 29.

P.S. Nobody thinks that this is a substitute for actually building something!
