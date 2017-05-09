package com.rkm;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;

import com.rkm.de.SolveByDiffEvo;

import arduino.LowerMoments;

public class ReferenceSimulation {
	static int seed = 123456;
	static RandomGenerator rand = new MersenneTwister(seed);
	public static int NumDiffs = 6;
	static double[] measure = new double[NumDiffs+2]; // TODO: bad practice, gets around DifferentialEvolution limitations
	
	double V1, V2;
	VoltageRef v1 = new VoltageReference("LTC6655B-5 1", (V1=5.00020), 1.66, 2);
	VoltageRef v2 = new VoltageReference("LTC6655B-5 2", (V2=4.99985), 1.66, 2);
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
		    	{ 1,  1,  0,  0}, // not measurements--must sum to ~10
		    	{ 0,  0,  1,  1}}, false);
	
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
	private double[] getInitialSum() {
		double[] sum = new double[2];
		sum[0] = V1 + V2;
		sum[1] = R2a+R2b + (R1a+R1b)/2;
		return sum;
	}	/**
	 * Simulate a measurement by first gathering voltages from simulated references.
	 * These voltage are not actually know to measurement circuit so be careful in
	 * exposing these values.
	 */	
	public void solve(double temperature, double[] prevSum) {
		//double[] measure = new double[NumDiffs+2];
		//
		// Ask each of the references where they are "at":
		//
		final int VALS = 3;
		double[] val1 = new double[VALS];
		v1.next(temperature, val1);
		double[] val2 = new double[VALS];
		v2.next(temperature, val2);
		double[] val3 = new double[VALS];
		max6350.next(temperature, val3);
		double[] val4 = new double[VALS];
		lt6654.next(temperature, val4);
		//
		// Send back appropriate deltas as measurements
		int i = 0;
		
		measure[i++] = val1[0]-val2[0];
		measure[i++] = val1[1]-val3[0];
		measure[i++] = val1[2]-val4[0];
		measure[i++] = val2[1]-val3[0];
		measure[i++] = val2[2]-val4[0];
		measure[i++] = val3[0]-val4[0];
		measure[i++] = prevSum[0]; // imputed from previous solution--maybe should be averaged (until temperature changes)
		measure[i++] = prevSum[1]; // imputed

		//solveUsingDiffEvo(); // Not satisfactory		
		
		RealVector data = new ArrayRealVector(measure, false); // false meaning not copied, should it be??
		RealVector solution = solver.solve(data);
		
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g\n",
				solution.getEntry(0),solution.getEntry(1),solution.getEntry(2),solution.getEntry(3));
		// Second solution
		/*solution = solver2.solve(data);
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g\n",
				solution.getEntry(0),solution.getEntry(1),solution.getEntry(2),solution.getEntry(3));*/
		System.out.printf("delta=%g \t dv2= %g \t dr3 = %g \t dr2= %g\n\n", 
				solution.getEntry(0)-mean(val1),
				solution.getEntry(1)-mean(val2),
				solution.getEntry(2)-mean(val3),
				solution.getEntry(3)-mean(val4));
		//
		// As promised, returning the current solution sums to be used in next iteration.
		//
		prevSum[0]= solution.getEntry(0) + solution.getEntry(1);
		prevSum[1]= solution.getEntry(2) + solution.getEntry(3);
		// Experiment: what happens if we use the exact values?? 
		// If this (could) be done the error is in the .1-.01 microvolt range. Rather than the
		// typical 3-20  microvolt range.
		// prevSum[0]= val1[0]+val2[2];
		// prevSum[1]= val3[1]+val4[0];
	}
	
	private void solveUsingDiffEvo() {
		SolveByDiffEvo solver = new SolveByDiffEvo();
		double[] solution = solver.getBest();
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g\n",
				solution[0], solution[1], solution[2], solution[3]);
		
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
		final int MAX = 150; // up 5, down 5
		double maxTemp = VoltageRef.RoomTemp;;
		for (int cycle = 0; cycle < 80; cycle++) {
			double temp = VoltageRef.RoomTemp;
			for (int i = 0; i < MAX; i++) {
				System.out.println("Temperature = "+temp);
				for (int k = 0; k < 10; k++) {
					/*double[] refs = sim.actual(temp);
					double[] m1 = sim.oneMeasurement(refs, sum);*/
					sim.solve(temp, sum);

				}
				temp += (i>=MAX/2)? -0.15 : 0.15;
				maxTemp = Math.max(temp, maxTemp);
				if (i==MAX/2) {
					System.out.println("max");
				}
			}
		}
		System.out.println("Maximum temperature = "+maxTemp);

	}
	
	/**
	 * Used to pass latest measurement to differential evolution classes
	 * @return measurement vector
	 */
	public static double[] getLastMeasurement() {
		return measure;
	}

}
