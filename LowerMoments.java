package arduino;
/*
* Copyright (c) 2003, the JUNG Project and the Regents of the University
* of California
* All rights reserved.
*
* This software is open-source under the BSD license.
* See http://jung.sourceforge.net/license.txt for a description.
*/

public class LowerMoments {
    /**
     * A data structure representing the (lower) central moments of a distribution including: <ul>
     * <li> the mean </li>
     * <li> the variance </li>

     * Data values are each passed into this data structure via the accumulate(...) method
     * and the corresponding central moments are updated on each call
     * 
     * Elided skewness and kurtosis. RKM 4/3/2017
     *
     * @author Didier H. Besset (modified by Scott White and, subsequently, by Leland Wilkinson)
     */
    private double[] moments;

    public LowerMoments() {
        moments = new double[3];
    }

    public void accumulate(double x) {
        if (Double.isNaN(x) || Double.isInfinite(x))
            return;
        double n = moments[0];
        double n1 = n + 1;
        //double n2 = n * n;
        double delta = (moments[1] - x) / n1;
        double d2 = delta * delta;
        //double d3 = delta * d2;
        double r1 = n / n1;
        //moments[4] += 4 * delta * moments[3] + 6 * d2 * moments[2] + (1 + n * n2) * d2 * d2;
        //moments[4] *= r1;
        //moments[3] += 3 * delta * moments[2] + (1 - n2) * d3;
        //moments[3] *= r1;
        moments[2] += (1 + n) * d2;
        moments[2] *= r1;
        moments[1] -= delta;
        moments[0] = n1;
    }

    public double mean() {
        return moments[1];
    }

    public double count() {
        return moments[0];
    }

    public double standardDeviation() {
        return Math.sqrt(variance());
    }

    public double variance() {
        if (moments[0] < 2)
            return Double.NaN;
        return moments[2] * moments[0] / (moments[0] - 1);
    }

	public void clear() {
		moments[0] = moments[1] = moments[2] = 0;		
	}
	
	@Override public String toString() {
	    StringBuilder result = new StringBuilder();
	    result.append(this.mean());
	    result.append(" sd=");
	    result.append(this.standardDeviation());
		return result.toString();	    
	}
}
