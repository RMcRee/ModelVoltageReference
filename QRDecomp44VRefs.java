/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package arduino;
import org.apache.commons.math3.util.FastMath;
/**
 * Stole this excerpt from QRDecomposition in commons.math. Specialized to solve only the
 * coefficient matrix for four (44, get it?) differential voltage references.
 * We can now bury this into a microcontroller w/out having all the baggage. 
 * @author randallm
 *
 */
public class QRDecomp44VRefs {

	double[] rDiag = {-2.0, -2.0, -1.8708286933869707, 1.8516401995451028};
	double[][] qrt = {
		{3.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0},
		{0.0, 2.3333333333333335, 0.3333333333333333, 1.0, 1.0, 0.0, 1.3333333333333333, 0.0}, 
		{0.5, 0.5000000000000003, 2.2279715505298276, -0.4285714285714285, 0.5714285714285715, 1.0, 0.9285714285714286, 1.0},
		{0.5, 0.5, 0.2672612419124243, -1.911332872020394, -0.587076436699612, -0.527383764224321, 0.7960007903631305, 1.472616235775679}};
    
	public double[] solve(double[] y) {
        final int n = qrt.length;
        final int m = qrt[0].length;

        final double[] x = new double[n];

        // apply Householder transforms to solve Q.y = b
        for (int minor = 0; minor < FastMath.min(m, n); minor++) {

            final double[] qrtMinor = qrt[minor];
            double dotProduct = 0;
            for (int row = minor; row < m; row++) {
                dotProduct += y[row] * qrtMinor[row];
            }
            dotProduct /= rDiag[minor] * qrtMinor[minor];

            for (int row = minor; row < m; row++) {
                y[row] += dotProduct * qrtMinor[row];
            }
        }

        // solve triangular system R.x = y
        for (int row = rDiag.length - 1; row >= 0; --row) {
            y[row] /= rDiag[row];
            final double yRow = y[row];
            final double[] qrtRow = qrt[row];
            x[row] = yRow;
            for (int i = 0; i < row; i++) {
                y[i] -= yRow * qrtRow[i];
            }
        }

        return x;
    }
	
}
