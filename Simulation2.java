package com.rkm;

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

public class Simulation2 {
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
	 * Keep average, stddev of all four references, v1, v2, r1, r2
	 */
	LowerMoments averageOf4 = new LowerMoments();
	LowerMoments[] solnAverage = new LowerMoments[]{new LowerMoments(),new LowerMoments(),
			                                        new LowerMoments(),new LowerMoments()};
	double maxDelta = 0;
	/**
	 * The coefficients for our least squares problem:
	 *  v1 + v2 +  r1 +  r2 = measure[0]  --this is needed or LSQ algo picks arbitrary values; too under-constrained
	 *  v1 - v2 + 0r1 + 0r2 = measure[1]
	 * ...
	 * 0v1+ 0v2 +  r1 - r2  = measure[6]
	 */
	
	// 10 solvers: 6 for all pairs of v1/v2...r1/r2 and v1 alone, v2,r1,r2 alone.
	private DecompositionSolver[] solver = new DecompositionSolver[10];
	//private DecompositionSolver solver2;
	
	public Simulation2() {
		RealMatrix coefficients = new Array2DRowRealMatrix(new double[][] {		    	
			{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
			{ 0,  1,  1,  1}}, true);
		solver[0] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  0,  1,  1}}, true);
		solver[1] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  1,  0,  1}}, true);
		solver[2] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  1,  1,  0}}, true);
		solver[3] = new QRDecomposition(coefficients).getSolver(); 
		// All pairs
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  1,  0,  0}}, true);
		solver[4] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  0,  1,  0}}, true);
		solver[5] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 1,  0,  0,  1}}, true);
		solver[6] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  1,  1,  0}}, true);
		solver[7] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  1,  0,  1}}, true);
		solver[8] = new QRDecomposition(coefficients).getSolver(); 
		coefficients = new Array2DRowRealMatrix(new double[][] {		    	
	    	{ 1, -1,  0,  0}, { 1,  0, -1,  0},{ 1,  0,  0, -1},{ 0,  1, -1,  0},{ 0,  1,  0, -1},{ 0,  0,  1, -1},
	    	{ 0,  0,  1,  1}}, true);
		solver[9] = new QRDecomposition(coefficients).getSolver(); 
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
		measure[i++] = val1[0]-val2[0]; // v1 - v2  0
		measure[i++] = val1[1]-val3[0]; // v1 - r1  1
		measure[i++] = val1[2]-val4[0]; // v1 - r2  2
		measure[i++] = val2[1]-val3[1]; // v2 - r1  3
		measure[i++] = val2[2]-val4[1]; // v2 - r2  4
		measure[i++] = val3[2]-val4[2];	// r1 - r2	5
		
		RealVector solution[] = new RealVector[solver.length]; 
		// brute force: generate all solutions for all equations
		LowerMoments[] mean = generateSolutions(prevSolution, solution, measure);
		printit(temperature,  mean, val1, val2, val3, val4);
		//
		// As promised, returning the current, best, solution to be used in next iteration.
		//

	}
	
	private LowerMoments[] generateSolutions(double[] prevSolution, RealVector[] solution, double[] measure) {
		int worst = -1; // illegal
		double maxResidual = Double.MIN_VALUE;
		for (int i = 0; i < solver.length; i++) {
			measure[6] = equationSelect(i, prevSolution);
			RealVector data = new ArrayRealVector(measure, true);
			solution[i] = solver[i].solve(data);
			if (i <= 3) {
				double residual = getResiduals(solution[i], prevSolution);
				if (maxResidual < residual) {
					maxResidual = residual;
					worst = i;
				}
			}
		}
		boolean[] include = excludeOneMember(worst, maxResidual);
		LowerMoments[] mean = new LowerMoments[prevSolution.length];
		for (int j = 0; j < mean.length; j++) { mean[j] = new LowerMoments(); }
		for (int i = 0; i < solver.length; i++) {
			if (include[i]) {
				// use solution[i] to contribute to prevsolution
				for (int j = 0; j < mean.length; j++) { 
					mean[j].accumulate(solution[i].getEntry(j)); 
				}
			}
		}
		//
		// As promised, returning the current, best, solution to be used in next iteration.
		//
		for (int j = 0; j < mean.length; j++) {
			prevSolution[j] = mean[j].mean(); 
		}
		accumulate(prevSolution, solnAverage);
		return mean;		
	}
	
	private double getResiduals(RealVector sol, double[] prevSolution) {
		double sum = 0;
		for (int i=0; i < prevSolution.length; i++) {
			sum += FastMath.abs(sol.getEntry(i) - prevSolution[i]);			
		}	
		return sum;
	}
	
	private double getResidualsSansOutlier(RealVector sol, double[] prevSolution) {		
		double max = Double.MIN_VALUE;
		int maxIndex = -1;
		for (int i=0; i < prevSolution.length; i++) {
			double delta = FastMath.abs(sol.getEntry(i) - prevSolution[i]);
			if (max < delta) {
				max = delta; 
				maxIndex = i; // need to know the index of the maximum deviation
			}
		}
		double sum = 0;
		for (int i=0; i < prevSolution.length; i++) {
			if (i != maxIndex) {
				sum += FastMath.abs(sol.getEntry(i) - prevSolution[i]);
			}
		}	
		return sum;
	}
	
	private double equationSelect(int i, double[] prevSolution) {
		switch (i) {
		case 0: return prevSolution[1]+prevSolution[2]+prevSolution[3];
		case 1: return prevSolution[0]+prevSolution[2]+prevSolution[3];
		case 2: return prevSolution[0]+prevSolution[1]+prevSolution[3];
		case 3: return prevSolution[0]+prevSolution[1]+prevSolution[2];
		case 4: return prevSolution[0] + prevSolution[1];
		case 5: return prevSolution[0] + prevSolution[2];
		case 6: return prevSolution[0] + prevSolution[3];
		case 7: return prevSolution[1] + prevSolution[2];
		case 8: return prevSolution[1] + prevSolution[3];
		case 9: return prevSolution[2] + prevSolution[3];
		}
		throw new IllegalArgumentException();
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
	
	private void printit(double temperature, LowerMoments[] mean, 
			             double[] val1, double[] val2, double[] val3, double[] val4) {
		double v1 = mean(val1), v2 = mean(val2), r1=mean(val3), r2 = mean(val4);
		double average = (v1+v2+r1+r2)/4;
		averageOf4.accumulate(average);
		//
		// Only print when things get worse, to save on eyesight
		//
		if (maxDelta < mean[0].mean()-v1) {
			maxDelta = mean[0].mean()-v1;
		} else {
			//return;
		}
		if (FastMath.abs(mean[0].mean()-v1)>1e-5){
			System.out.println("oops");
		}
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g \tavg=%10.7g\t<-lsq soln\n",
				mean[0].mean(),mean[1].mean(),mean[2].mean(),mean[3].mean(), average);
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g \t\t\t<-actual\n",
				v1,v2,r1,r2);
		// Second solution
		/*solution = solver2.solve(data);
		System.out.printf("v1  = %10.7g \t v2 = %10.7g \t r1  = %10.7g \t r2 = %10.7g\n",
				solution.getEntry(0),solution.getEntry(1),solution.getEntry(2),solution.getEntry(3));*/
		// Calculate average of all four
		System.out.printf("delta=%g \t dv2= %g \t dr3 = %g \t dr2= %g\t temperature=%3.2g\n\n", 
				mean[0].mean()-v1,
				mean[1].mean()-v2,
				mean[2].mean()-r1,
				mean[3].mean()-r2, temperature);		
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
		Simulation2 sim = new Simulation2();
		
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
