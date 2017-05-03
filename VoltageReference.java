package com.rkm;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

import arduino.LowerMoments;
import arduino.PinkNoise;

/**
 * Simulate a voltage reference.
 * It has a nominal voltage
 * It has a temperature coefficient, expressed in ppm per degrees kelvin
 * It has a pink noise source, expressed in ?
 * We do *not* model white noise.
 * We do *not* model long term drift
 * 
 * Time scale is assumed to be one-tenth of a second. E.g. each call to next() is assumed to be
 * one-tenth of a second after the previous one. This only matters if long-term drift is modeled.
 * 
 * @author Randall McRee
 * 4/30/2017
 *
 */
public class VoltageReference {

	private static final double Alpha = 1.1;
	private static final int Poles = 10;
	
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
	 * Simple statistics: mean and standard deviation
	 */
	protected LowerMoments stats;
	/**
	 * Generator for random numbers. Also used by embedded pink noise object.
	 */
	protected RandomGenerator rndom;
	
	public VoltageReference(String _whoami, double value) {
		this(_whoami, value, 5*ppm.one, 5, 0l);
	}
	/**
	 * Create voltage reference with given (statistical) properties
	 * @param _whoami  Identity to use when printing
	 * @param nominal  Nominal voltage
	 * @param noise_pp Flicker noise
	 * @param _tempco  Temperature coefficient. PPM per degrees centigrade.
	 * @param seed
	 */
	public VoltageReference(String _whoami, double nominal, double noise_pp, double _tempco, long seed) {
		this.whoami = _whoami;
		rndom = (seed==0)? new MersenneTwister() : new MersenneTwister(seed);
		noise = new PinkNoise(Alpha, Poles, rndom);
		noise_scale = noise_pp * ppm.one; // TODO: correct this constant
		stats = new LowerMoments();
	}
	
	public double next(double temperature) {
		double value = voltage + ( noise.nextValue()*noise_scale );
		stats.accumulate(value);
		return value;
	}
	
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
}
