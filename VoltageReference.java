package com.rkm;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

import arduino.LowerMoments;
import arduino.PinkNoise;

/**
 * Simulate a voltage reference.
 * It has a nominal voltage.
 * It has a temperature coefficient, expressed in ppm per degrees kelvin.
 * It has a pink noise source, expressed in micro-volts peak-to-peak (prevalent in datasheets).
 * We do *not* model white noise.
 * We do *not* model long term drift.
 * 
 * Time scale is ambiguous.
 * 
 * @author Randall McRee
 * 4/30/2017
 *
 */
public class VoltageReference implements VoltageRef {

	private static final double Alpha = 1.0;
	private static final int Poles = 10;
	private static final double pp_to_rms_1_0 = 9.3080691;
	private static final double pp_to_rms_1_1 = 10.2289;
	private static final double DriftThreshold = 0.999; // 1.0 is no drift. 0.998 seems to be realistic (depending on your timescale).
	
	
	/**
	 * The identity of the reference. Used in toString(), only.
	 */
	protected String whoami;
	/**
	 * Nominal voltage. Modified by time, temperature and noise source.
	 */
	protected double voltage;
	/**
	 * Pink (flicker) noise of the correct magnitude.
	 */
	protected PinkNoise noise;
	protected double noise_scale; 
	/**
	 * Temperature coefficient in ppm-per-degree centigrade
	 */
	protected double tempco;
	/**
	 * Our internal temperature. It lags ambient temperature.
	 * RKM made this up, it is supposed to model the fact that references are
	 * obviously sensitive to temperature and each has its own unique environment.
	 */
	protected double internalTemp;
	/**
	 * Simple statistics: mean and standard deviation
	 */
	protected LowerMoments stats;
	/**
	 * Generator for random numbers. Also used by embedded pink noise object.
	 */
	protected RandomGenerator rndom;
	
	/**
	 * Ctor for subclasses, averaging and stacked type references.
	 */
	protected VoltageReference() {
		
	}
	
	public VoltageReference(String _whoami, double value) {
		this(_whoami, value, 5*ppm.one, 5, 0l);
	}
	/**
	 * Create voltage reference with given (statistical) properties. Random generator seed is
	 * different for each instantiation.
	 * 
	 * @param _whoami  Identity to use when printing
	 * @param nominal  Nominal voltage
	 * @param noise_pp Flicker noise in units of microvolts peak-peak
	 * @param _tempco  Temperature coefficient. PPM per degrees centigrade.
	 */	
	public VoltageReference(String _whoami, double nominal, double noise_pp, double _tempco) {
		this(_whoami, nominal, noise_pp, _tempco, 0);
	}
	/**
	 * Create voltage reference with given (statistical) properties
	 * @param _whoami  Identity to use when printing
	 * @param nominal  Nominal voltage
	 * @param noise_pp Flicker noise in units of microvolts peak-peak
	 * @param _tempco  Temperature coefficient. PPM per degrees centigrade.
	 * @param seed
	 */
	public VoltageReference(String _whoami, double nominal, double noise_pp, double _tempco, long seed) {
		this.whoami = _whoami;
		voltage = nominal;
		rndom = (seed==0)? new MersenneTwister() : new MersenneTwister(seed);
		noise = new PinkNoise(Alpha, Poles, rndom);
		noise_scale = (noise_pp  * ppm.one)/(pp_to_rms_1_0);
		tempco = _tempco * voltage * ppm.one;
		internalTemp = RoomTemp;
		stats = new LowerMoments();
	}
	/* (non-Javadoc)
	 * @see com.rkm.VoltageRef#next(double)
	 */
	@Override
	public double next(double temperature) {
		drift();
		double value = voltage + ( noise.nextValue()*noise_scale);
		double deltaT = (temperature - internalTemp) * rndom.nextDouble(); // see discussion above
		internalTemp += deltaT;
		double deltaV = tempco * (internalTemp - RoomTemp);
		value += deltaV; // Temperature effect. Straight line.
		stats.accumulate(value);
		return value;
	}
	/* (non-Javadoc)
	 * @see com.rkm.VoltageRef#next(double, double[])
	 */
	@Override
	public void next(double temperature, double[] values) {
		drift();
		double deltaT = (temperature - internalTemp) * rndom.nextDouble(); // see discussion above
		internalTemp += deltaT;
		double deltaV = tempco * (internalTemp - RoomTemp);
		for (int i = 0; i < values.length; i++) {
			double value = voltage + ( noise.nextValue()*noise_scale);
			value += deltaV; // Temperature effect. Straight line. 
			stats.accumulate(value);
			values[i] = value;
		}		
	}
	
	/* (non-Javadoc)
	 * @see com.rkm.VoltageRef#toString()
	 */
	@Override public String toString() {
	    StringBuilder sb = new StringBuilder();
	    sb.append(whoami);
	    sb.append("=");
	    sb.append(voltage);
	    sb.append("@");
	    sb.append(tempco);
	    sb.append("ppm-per-K");
	    return sb.toString();
	}
	
	/* (non-Javadoc)
	 * @see com.rkm.VoltageRef#sd()
	 */
	@Override
	public double sd() {
		return stats.standardDeviation();
	}
	
	/**
	 * Experimental. Popcorn noise causes random shifts in output voltage.
	 * Currently not used.
	 * @return
	 */
	boolean popcorn() {
		if (rndom.nextDouble()>0.995) {
			return true;
		} else {
			return false;
		}
	}
	
	void drift() {
		if (rndom.nextDouble()>DriftThreshold){
			voltage += (voltage*ppm.one) * (rndom.nextDouble()-0.50);
		} 		
	}

}
