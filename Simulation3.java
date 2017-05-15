package com.rkm;


import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

import arduino.LowerMoments;

public class Simulation3 {
	private static final double StdDevThreshold = 1.7e-7;

	static long seed[] = new long[]{ 295075147L, 141650939L, 715225739L, 961748927L};
            
	static RandomGenerator rand = new MersenneTwister(414941729579l);
	public static int NumDiffs = 6;
	
	double V1, V2;
	VoltageRef lt6655_1 = new VoltageReference("LTC6655B-5 1", (V1=5.00020), 1.66, 2/*, seed[0]*/);
	VoltageRef lt6655_2 = new VoltageReference("LTC6655B-5 2", (V2=4.99960), 1.66, 2/*, seed[1]*/);
	// 
	// Do not use averaged or stacked references...obscures effects of drift.
	//
	double R1;
	VoltageReference refC = new VoltageReference("MAX6350 C", (R1=5.00043), 3, 1/*,seed[2]*/);
	//VoltageReference refD = new VoltageReference("MAX6350 D", (R1b=5.00022), 3, 1);
	VoltageRef max6350 = refC; // new AverageReference(refC, refD);
	
	double R2;
	VoltageReference refE = new VoltageReference("LT6654 E", (R2=5.00015), 3.7, -3/*,seed[3]*/);
	//VoltageReference refF = new VoltageReference("LT6654 F", (R2b=2.50018), 3.7, -3);
	VoltageRef lt6654 = refE; // new StackedReference(refE, refF);
	LowerMoments meanSum = new LowerMoments();
	/**
	 * Keep averages, stddev of all four references, v1, v2, r1, r2; also v1+v2+r1+r2
	 */
	LowerMoments   averageOf4  = new LowerMoments();
	LowerMoments[] solnAverage = new LowerMoments[]{new LowerMoments(),new LowerMoments(),
			                                        new LowerMoments(),new LowerMoments()};
	double maxDelta = 0;
	
	//  solvers: 
	private DecompositionSolver[] solver = new DecompositionSolver[8];
	//private DecompositionSolver solver2;
	
	public Simulation3() {
		// solve for three vrefs, x-y, x-z, y-z, x+y+z
		RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] {		    	
			{ 1, -1,  0}, { 1,  0, -1},{ 0,  1, -1},{ 1,  1,  1}}, true);
		solver[0] = new QRDecomposition(coefficients).getSolver();
		// Equally weighted (no drift) matrix
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0}, { 1,  0,  0, -1}, { 0,  1, -1,  0},
	    	{ 0,  1,  0, -1}, { 0,  0,  1, -1}, { 1,  1,  1,  1}}, true);
		solver[1] = new QRDecomposition(coefficients).getSolver(); 
		//
		// All pairs, pick the best pair to use...
		//
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  1,  0,  0}}, true);
		solver[2] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  0,  1,  0}}, true);
		solver[3] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  0,  0,  1}}, true);
		solver[4] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  1,  1,  0}}, true);
		solver[5] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  1,  0,  1}}, true);
		solver[6] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  0,  1,  1}}, true);
		solver[7] = new QRDecomposition(coefficients).getSolver(); 		
	}
	/**
	 * Initialize previous solution array. This 'seeds' the initial solution.
	 * @return array maintained by solver consisting of current, imputed values.
	 */
	double[] getInitialValues() {
		double[] soln = new double[4];
		soln[0] = V1;
		soln[1] = V2;
		soln[2] = R1;
		soln[3] = R2;
		accumulate(soln, solnAverage);
		return soln;
	}
	
	private void accumulate(double[] soln, LowerMoments[] average) {
		for (int i = 0; i < soln.length; i++) {
			average[i].accumulate(soln[i]);
		}		
	}
	/**
	 * Simulate a measurement by first gathering voltages from simulated references.
	 * These voltage are not actually known to measurement circuit so be careful in
	 * exposing these values.
	 */	
	public void solve(double temperature, double[] prevSolution) {
		//
		// Ask each of the references where they are "at":
		//
		final int VALS = 3;
		double[] val1 = new double[VALS];
		lt6655_1.next(temperature, val1); // Get three measurements at once, for each reference in turn. 
		// Experiment: make this one odd man out
		//val1[0]+= 5e-6; val1[1]+= 5e-6; val1[2]+= 5e-6;
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
		
		// Find two references that have drifted least and use that pair to solve
		RealVector soln = generateSolutions(prevSolution, measure);
		printit(temperature, soln, val1, val2, val3, val4);
	}
	
	private RealVector generateSolutions(double[] prevSolution, double[] measure) {
		int worst = -1; // illegal
		double maxResidual = Double.MIN_VALUE;
		LowerMoments[] mean3 = new LowerMoments[prevSolution.length];
		for (int j = 0; j < mean3.length; j++) { mean3[j] = new LowerMoments(); }
		RealVector[] solution = new RealVector[4];
		for (int i = 0; i < 4; i++) {			
			RealVector data = equationSelect(i, measure, prevSolution);		
			solution [i] = solver[0].solve(data);
			switch (i) {
			case 0: //v1v2r1
				mean3[0].accumulate( solution[i].getEntry(0));
				mean3[1].accumulate( solution[i].getEntry(1));
				mean3[2].accumulate( solution[i].getEntry(2));
				break;
			case 1: // v1v2r2
				mean3[0].accumulate( solution[i].getEntry(0));
				mean3[1].accumulate( solution[i].getEntry(1));
				mean3[3].accumulate( solution[i].getEntry(2));
				break;
			case 2: // v1r1r2
				mean3[0].accumulate( solution[i].getEntry(0));
				mean3[2].accumulate( solution[i].getEntry(1));
				mean3[3].accumulate( solution[i].getEntry(2));				
				break;
			case 3: // v2r1r2
				mean3[1].accumulate( solution[i].getEntry(0));
				mean3[2].accumulate( solution[i].getEntry(1));
				mean3[3].accumulate( solution[i].getEntry(2));
				break;
			}		
		}
		int index = findSolution(mean3, solution, prevSolution);
		measure[6] = sum(prevSolution, index);
		RealVector data = new ArrayRealVector(measure, true);
		RealVector soln = solver[index].solve(data);					
		//
		// As promised, returning the current, best, solution to be used in next iteration.
		//
		for (int j = 0; j < prevSolution.length; j++) {
			prevSolution[j] = soln.getEntry(j); 
		}
		return soln;
	}
	private double sum(double[] prevSolution, int index) {
		switch(index) {
		case 1: return prevSolution[0]+prevSolution[1]+prevSolution[2]+prevSolution[3];
		case 2:	return prevSolution[0]+prevSolution[1];
		case 3:	return prevSolution[0]+prevSolution[2];
		case 4:	return prevSolution[0]+prevSolution[3];
		case 5:	return prevSolution[1]+prevSolution[2];
		case 6:	return prevSolution[1]+prevSolution[3];
		case 7:	return prevSolution[2]+prevSolution[3];
		default: throw new IllegalArgumentException();
		}
	}
	/**
	 * Helper method: return array sum
	 * @param array
	 * @return  Sum of all elements;
	 */
	private double sum(double[] array) {
		double sum = 0;
		for (int i = 0; i < array.length; i++) sum += array[i];
		return sum;
	}
	/**
	 * Return index of the solver that should be used to fully solve this iteration. Two of the references are
	 * most trusted and contribute fully, therefore, to final solution. Other two are inferred from measurement.
	 * @param mean
	 * @param solution
	 * @param prevSolution
	 * @return
	 */
	private int findSolution(LowerMoments[] mean, RealVector[] solution, double[] prevSolution) {
		int index = 0;
		for (int i = 0; i < mean.length; i++) {
			mean[i].accumulate(prevSolution[i]);
		}
		TreeMap<Double,Integer> pair = new TreeMap<Double,Integer>();
		for (int i = 0; i < mean.length; i++) {
			pair.put(mean[i].standardDeviation(), i);
		}
		if (pair.size() != 4) throw new IllegalArgumentException(); // could happen if two stddevs are exacty equal!
		int lowest = pair.firstEntry().getValue();
		double val = pair.firstKey();
		pair.remove(val);
		if (FastMath.abs(val-pair.lastKey()) < StdDevThreshold) {
			return 1; // use all sums solution
		}
		int low = pair.firstEntry().getValue();
		switch (lowest) {
		case 0: // v1
			switch (low) {
			case 0: throw new IllegalArgumentException(); 
			case 1: return 2; // v1v2
			case 2: return 3; // v1r1
			case 3: return 4; // v1r2
			}
		case 1: // v2
			switch (low) {
			case 0: return 2; // v1v2
			case 1: throw new IllegalArgumentException();
			case 2: return 5; // v2r1
			case 3: return 6; // v2r2
			}
		case 2: // r1
			switch (low) {
			case 0: return 3; // v1r1
			case 1: return 5; // v2r1
			case 2: throw new IllegalArgumentException();
			case 3: return 7; // r1r2
			}
		case 3: // r2
			switch (low) {
			case 0: return 4; // v1r2
			case 1: return 6; // v2r2
			case 2: return 7; // r1r2
			case 3: throw new IllegalArgumentException();
			}
		}
		throw new IllegalArgumentException();
	}
	// TODO: now we have permutation of three variables, not four each time
	private double getResiduals(RealVector sol, double[] prevSolution) {
		double sum = 0;
		for (int i=0; i < prevSolution.length; i++) {
			sum += FastMath.abs(sol.getEntry(i) - prevSolution[i]);			
		}	
		return sum;
	}
		
	private RealVector equationSelect(int i, double[] measure, double[]  prevSolution) {
		double[] meas = new double[4];
		switch (i) {
		case 0: //v1v2r1
			meas[0]= measure[0]; //v1-v2
			meas[1]= measure[1]; //v1-r1
			meas[2]= measure[3]; //v2-r1
			meas[3]= prevSolution[0]+prevSolution[1]+prevSolution[2];
			break;
		case 1: // v1v2r2
			meas[0]= measure[0]; //v1-v2
			meas[1]= measure[2]; //v1-r2
			meas[2]= measure[4]; //v2-r2
			meas[3]= prevSolution[0]+prevSolution[1]+prevSolution[3];
			break ;
		case 2: // v1r1r2
			meas[0]= measure[1]; //v1-r1
			meas[1]= measure[2]; //v1-r2
			meas[2]= measure[5]; //r1-r2
			meas[3]= prevSolution[0]+prevSolution[2]+prevSolution[3];
			break ;
		case 3: // v2r1r2
			meas[0]= measure[3]; //v2-r1
			meas[1]= measure[4]; //v2-r2
			meas[2]= measure[5]; //r1-r2
			meas[3]= prevSolution[1]+prevSolution[2]+prevSolution[3];
			break;

		}
		return new ArrayRealVector(meas, false);
	}
	
	private boolean[] excludeOneMember(int which, double maxResidual) {
		if (maxResidual < 5e-6)
			return new boolean[] {true, true, true, true, true, true, true, true,  true,  true};
		switch (which) {
		//                             v1     v2    r1    r2    v1v2   v1r1   v1r2   v2r1   v2r2   r1r2
		case 0: return new boolean[] {true, false, false, false, false, false, false, true,  true,  true};  // exclude v1
		case 1: return new boolean[] {false, true, false, false, false, true,  true,  false, false, true};  // exclude v2
		case 2: return new boolean[] {false, true, true, false, false,  false, true,  false, true,  false}; // exclude r1
		case 3: return new boolean[] {false, false, false, true, true,  true,  false, true,  false, false}; // exclude r2
		}
		throw new IllegalArgumentException();
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
			//return;
		}
		if (FastMath.abs(solution.getEntry(0)-v1)>=1e-5) {
			System.out.println("oh");
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
		Simulation3 sim = new Simulation3();
		
		System.out.println("Simple modeling of four voltage references and their differential measurements");
		double[] sum = sim.getInitialValues();		 
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
				for (int k = 0; k < 10; k++) {
					/*double[] refs = sim.actual(temp);
					double[] m1 = sim.oneMeasurement(refs, sum);*/
					sim.solve(temp, sum);
				}
				
				temp += (i>=MAX/2)? -0.05 : 0.05; // ramp up then down, over and over 				
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
