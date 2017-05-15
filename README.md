# ModelVoltageReference
Model pink noise and temperature dependence of voltage references

Short discussion of my purpose in creating this here:
http://www.eevblog.com/forum/metrology/software-simulation-of-a-voltage-reference-noise-tempco-long-term-drift/

The idea is to simulate some of the more salient properties of voltage references that affect their performance,
namely noise and temperature.

For example, let's say that you build a circuit which can measure the difference between two voltage references and 
amplify that difference (e.g. an instrumentation amplifier). The simulation ought to give us some insight into the
noise and temperature dependence of the INAMP output. Furthermore, if we have a set of, say, four references, and know their starting values can we track their changes over temperature and long-term drift using just differential measurements?

Simulation shows that differential-only measurements may be a viable way to accurately track (not measure) a set of references. After developing this, I googled "NIST least squares metrology" and came up with this 1974 paper which discusses measuring differences in gage blocks (coming up with an almost identical least squares equation):
http://emtoolbox.nist.gov/Publications/NBSIR74-587.asp
See especially appendix at end. Note that author talks about adding constraint to LSQ: (A+B)=K on page 29.

You do need to start with a somewhat accurate measurement (this is initialSum in the simulation). The tracking accuracy will always be with respect to this initial accuracy. You do not get something for nothing. I think this is still valuable since readily available voltage references differ a fair amount with temperature and this scheme should allow us to quantify that.

After adding the possibility of drift it seems that the method is not a panacea. Although there may be some interesting applications, it seems that realistic drift and temperature rates cause the lsq solution to deviate from reality by an arbitrary amount.

Simulation3 does show a way in which a single reference can be detected to be erring and, if the other references have not drifted, it can be corrected (meaning the lsq solution can be updated to have the correct reading). This may be interesting by itself.

