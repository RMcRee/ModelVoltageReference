package com.rkm;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

import arduino.LowerMoments;

public class ReferenceSimulation {
	static int seed = 123456;
	static RandomGenerator rand = new MersenneTwister(seed);
	
	public static void main(String[] args)
	{
		System.out.println("Simple modeling of voltage references");
		double temp = 23; // degrees centigrade
		LowerMoments mom = new LowerMoments();
		VoltageReference ref1 = new VoltageReference("LTC6655-2.5", 2.5001, 0.83, 2);
		VoltageReference ref2 = new VoltageReference("LTC6655-2.5", 2.4998, 0.83, 2);
		for (int k = 0; k<20; k++) {			
			for (int i = 0; i < 1000000; i++) {
				double v = ref1.next(temp);
				double w = ref2.next(temp);
				double avg = (v+w)/2;
				mom.accumulate(avg);			
			}
			// sd stands for standard deviation. In this case, it should equate directly to rms noise.
			System.out.printf("%g \t%g \t %s\n", ref1.sd(), ref2.sd(), mom.toString());
			temp += 0.5;
		}
	}

}
