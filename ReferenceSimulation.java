package com.rkm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

import arduino.LowerMoments;

public class ReferenceSimulation {
	static int seed = 123456;
	static RandomGenerator rand = new MersenneTwister(seed);
	public static int NumDiffs = 6;
	
	double V1, V2;
	VoltageRef lt6655_1 = new VoltageReference("LTC6655B-5 1", (V1=5.00020), 1.66, 2);
	VoltageRef lt6655_2 = new VoltageReference("LTC6655B-5 2", (V2=4.99960), 1.66, 2);
	// Next two are averaged together
	double R1a, R1b;
	VoltageReference refC = new VoltageReference("MAX6350 C", (R1a=5.00061), 3, 1);
	VoltageReference refD = new VoltageReference("MAX6350 D", (R1b=5.00022), 3, 1);
	VoltageRef max6350 = new AverageReference(refC, refD);
	
	// Two stacked lt6654
	double R2a, R2b;
	VoltageReference refE = new VoltageReference("LT6654 E", (R2a=2.49997), 3.7, -3);
	VoltageReference refF = new VoltageReference("LT6654 F", (R2b=2.50018), 3.7, -3);
	VoltageRef lt6654 = new StackedReference(refE, refF);
	LowerMoments meanSum = new LowerMoments();
	/**
	 * Keep average, stddev of all four references, v1, v2, r1, r2
	 */
	LowerMoments averageOf4 = new LowerMoments();
	double maxDelta = Double.MIN_VALUE;
	/**
	 * The coefficients for our least squares problem:
	 *  v1 + v2 +  r1 +  r2 = measure[0]  --this is needed or LSQ algo picks arbitrary values; too under-constrained
	 *  v1 - v2 + 0r1 + 0r2 = measure[1]
	 * ...
	 * 0v1+ 0v2 +  r1 - r2  = measure[6]
	 */
	RealMatrix coefficients =
		    new Array2DRowRealMatrix(new double[][] {
		    	
		    	{ 1, -1,  0,  0}, 
		    	{ 1,  0, -1,  0},
		    	{ 1,  0,  0, -1},
		    	{ 0,  1, -1,  0},
		    	{ 0,  1,  0, -1},
		    	{ 0,  0,  1, -1},
		    	{ 1,  1,  0,  0} // not measurements--must sum to ~10
		    	}, false);
	
	private DecompositionSolver solver;
	//private DecompositionSolver solver2;
	
	public ReferenceSimulation() {
		solver = new QRDecomposition(coefficients).getSolver(); // faster than SVD
		//solver2 = new SingularValueDecomposition(coefficients).getSolver();		
	}
	/**
	 * Initialize sum array. This 'seeds' the initial solution.
	 * @return array maintained by solver consisting of current, imputed values.
	 */
	double[] getInitialSum() {
		double[] sum = new double[2];
		sum[0] = V1 + V2;
		sum[1] = R2a+R2b + (R1a+R1b)/2;
		return sum;
	}	/**
	 * Simulate a measurement by first gathering voltages from simulated references.
	 * These voltage are not actually known to measurement circuit so be careful in
	 * exposing these values.
	 */	
	public void solve(double temperature, double[] prevSum) {
		//
		// Ask each of the references where they are "at":
		//
		final int VALS = 3;
		double[] val1 = new double[VALS];
		lt6655_1.next(temperature, val1); // Get three measurements at once, for each reference in turn. 
		double[] val2 = new double[VALS];
		lt6655_2.next(temperature, val2); // Would be unrealistic to have one measurement stand in for 
		double[] val3 = new double[VALS];
		max6350.next(temperature, val3);  // three real measurements, since each has its own noise.
		double[] val4 = new double[VALS];
		lt6654.next(temperature, val4);   //
		//
		// Now calculate the actual measurements as the delta of the references
		int i = 0;
		double[] measure = new double[NumDiffs+1];
		measure[i++] = val1[0]-val2[0];
		measure[i++] = val1[1]-val3[0];
		measure[i++] = val1[2]-val4[0];
		measure[i++] = val2[1]-val3[1];
		measure[i++] = val2[2]-val4[1];
		measure[i++] = val3[2]-val4[2];
		measure[i++] = prevSum[0]; // imputed from previous solution
		//measure[i++] = prevSum[1]; // imputed	
		
		RealVector data = new ArrayRealVector(measure, false); // false meaning not copied, should it be??
		RealVector solution = solver.solve(data);
		
		printit(temperature, solution, val1, val2, val3, val4);
		//
		// As promised, returning the current solution sums to be used in next iteration.
		//
		prevSum[0]= solution.getEntry(0) + solution.getEntry(1);
		//prevSum[1]= solution.getEntry(2) + solution.getEntry(3);
	}
	
	private void printit(double temperature, RealVector solution, 
			             double[] val1, double[] val2, double[] val3, double[] val4) {
		double v1 = mean(val1), v2 = mean(val2), r1=mean(val3), r2 = mean(val4);
		double average = (v1+v2+r1+r2)/4;
		averageOf4.accumulate(average);
		//
		// Only print when things get worse, to save on eyesight
		//
		if (maxDelta < solution.getEntry(0)-v1) {
			maxDelta = solution.getEntry(0)-v1;
		} else {
			return;
		}
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g \tavg=%10.7g\t<-lsq soln\n",
				solution.getEntry(0),solution.getEntry(1),solution.getEntry(2),solution.getEntry(3), average);
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g \t\t\t<-actual\n",
				v1,v2,r1,r2);
		// Second solution
		/*solution = solver2.solve(data);
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g\n",
				solution.getEntry(0),solution.getEntry(1),solution.getEntry(2),solution.getEntry(3));*/
		// Calculate average of all four
		System.out.printf("delta=%g \t dv2= %g \t dr3 = %g \t dr2= %g\t temperature=%3.2g\n\n", 
				solution.getEntry(0)-v1,
				solution.getEntry(1)-v2,
				solution.getEntry(2)-r1,
				solution.getEntry(3)-r2, temperature);		
	}
	
	private double mean(double[] val) {
		double mean = 0;
		for (int i = 0; i < val.length; i++) {
			mean += val[i];
		}
		return mean/val.length;
	}	
	
	public static void main(String[] args)
	{		
		ReferenceSimulation sim = new ReferenceSimulation();
		
		System.out.println("Simple modeling of four voltage references and their differential measurements");
		double[] sum = sim.getInitialSum();		 
		final int MAX = 300; 
		double maxTemp = VoltageRef.RoomTemp;
		/*sim.solve( VoltageRef.RoomTemp, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);
		sim.solve( VoltageRef.RoomTemp+5, sum);*/

		for (int cycle = 0; cycle < 365; cycle++) {
			double temp = VoltageRef.RoomTemp;
			for (int i = 0; i < MAX; i++) {
				//System.out.println("Temperature = "+temp);
				maxTemp = Math.max(temp, maxTemp);
				for (int k = 0; k < 1000; k++) {
					/*double[] refs = sim.actual(temp);
					double[] m1 = sim.oneMeasurement(refs, sum);*/
					sim.solve(temp, sum);
				}
				
				temp += (i>=MAX/2)? -0.10 : 0.10; // ramp up then down, over and over 				
				if (i==MAX/2) {
					//System.out.println("max");
				}
			}
		}
		System.out.println("Maximum temperature = "+maxTemp);
		System.out.printf("Average-of-4: mean = %8.7g,  max = %8.7g, min=%8.7g, stdev=%g",
				sim.averageOf4.mean(), sim.averageOf4.getMax(), sim.averageOf4.getMin(), sim.averageOf4.standardDeviation());

	}

}
