package com.rkm;

import arduino.LowerMoments;

public class AverageReference implements VoltageRef {

	VoltageReference ref1;
	VoltageReference ref2;
	LowerMoments avg;
	
	public AverageReference(VoltageReference _r1, VoltageReference _r2) {
		ref1 = _r1;
		ref2 = _r2;
		avg = new LowerMoments();
	}

	@Override
	public double next(double temperature) {
		double v1 = ref1.next(temperature);
		double v2 = ref2.next(temperature);
		double av = (v1+v2)/2;
		avg.accumulate(av);
		return av;
	}

	@Override
	public void next(double temperature, double[] values) {
		double[] moreVals = new double[values.length];
		ref1.next(temperature, moreVals);
		ref1.next(temperature, values);
		for (int i = 0; i < values.length; i++) {
			values[i] = (moreVals[i] + values[i])/2;
		}
	}
	
	@Override
	public double sd() {
		return avg.standardDeviation();
	}
	
	@Override 
	public String toString() {
	    StringBuilder sb = new StringBuilder();
	    sb.append(ref1.whoami);
	    sb.append("+");
	    sb.append(ref2.whoami);
	    sb.append("=");
	    sb.append((ref1.voltage+ref2.voltage)/2);
	    sb.append("@");
	    sb.append(ref1.tempco+ref2.tempco);
	    sb.append("ppm-per-K");
	    return sb.toString();
	}

}
