package com.rkm;

public interface VoltageRef {
	/**
	 * Room temperature constant
	 */
	public static final double RoomTemp = 23; // degrees centigrade
	/**
	 * Return the voltage for the reference given the specified noise source and temperature.
	 * Of course, this is an approximation, no real references will ever behave exactly like this.
	 * 
	 * When the temperature changes it is assumed to only affect the reference by a random amount.
	 * This models the fact that all references, even though they have similar properties, also have
	 * thermal mass (and lag the ambient temp) as well as their own unique physical environment. So, we
	 * arbitrarily model this as a random change in temperature (wave hands). The change, although random,
	 * is always in the same direction as the given parameter.
	 * 
	 * @param temperature  Temperature in degrees centigrade. 23 degrees is the arbitrary starting point.
	 * @return voltage
	 */
	double next(double temperature);
	/**
	 * Like next(temp) but fills 'values' with measurements. All have noise. 
	 * Temperature does not change during measurement.
	 * 
	 * @param temperature  Temperature in degrees centigrade. 23 degrees is the arbitrary starting point.
	 * @param values       Array in which to return successive values from this reference
	 */
	void next(double temperature, double[] values);
	
	String toString();
	
	/**
	 * Standard deviation of measurements returned so far.
	 * @return
	 */
	double sd();

}