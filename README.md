# ModelVoltageReference
Model pink noise and temperature dependence of voltage references

Short discussion of my purpose in creating this here:
http://www.eevblog.com/forum/metrology/software-simulation-of-a-voltage-reference-noise-tempco-long-term-drift/

The idea is to simulate some of the more salient properties of voltage references that affect their performance,
namely noise and temperature.

For example, let's say that you build a circuit which can measure the difference between two voltage references and 
amplify that difference (e.g. an instrumentation amplifier). The simulation ought to give us some insight into the
noise and temperature dependence of the INAMP output.

Nobody thinks that this is a substitute for actually building something!
