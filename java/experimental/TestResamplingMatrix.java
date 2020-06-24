/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package experimental;

import core.ResamplingLens;
import java.util.ArrayList;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class TestResamplingMatrix {
    
    Random rng = new Random();
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        TestResamplingMatrix go = new TestResamplingMatrix();
    }

    public TestResamplingMatrix() {
        /**
         * Create a baseline data set and populate with random numbers.
         */
        Jama.Matrix baseline = new Jama.Matrix(10000, 4);
        for (int i = 0; i < baseline.getRowDimension(); i++) {
            for (int j = 0; j < baseline.getColumnDimension(); j++) {
                baseline.set(i, j, rng.nextGaussian());
            }
        } 

        boolean testMemory = true;
        if (testMemory) {
            testMemory(baseline);
        }

        boolean testResampling = true;
        if (testResampling) {
            testResampling(baseline);
        }
    }

    private void testResampling(Jama.Matrix baseline) {
        /**
         * Show that resampling works
         */
        System.out.println("Original");
        pmUtility.prettyPrint(baseline.getMatrix(0,9,0,3));
        System.out.println("Resampled");
        ResamplingLens lens = new ResamplingLens(baseline, rng.nextLong());
        pmUtility.prettyPrint(lens.getMatrix(0,9,0,3));
        
        System.out.println("Baseline submatrix:");
        pmUtility.prettyPrint(baseline.getMatrix(0, 0, 0, baseline.getColumnDimension()-1));
        System.out.println("Resampled submatrix:");
        pmUtility.prettyPrint(lens.getMatrix(0, 0, 0, 3));
    }

    private void testMemory(Jama.Matrix baseline) {
        /**
         * Compute memory usage using resampled matrices
         */
        System.out.println("ResampledMatrix");
        for (int numResampledMatrices = 100; numResampledMatrices <= 100 * 100; numResampledMatrices *= 10) {
            ArrayList<ResamplingLens> resampledList = new ArrayList<>();

            /**
             * Now, create a resampled data set using that baseline data as the
             * underlying foundation.
             */
            for (int r = 0; r < numResampledMatrices; r++) {
                resampledList.add(new ResamplingLens(baseline, rng.nextLong()));
            }
            System.out.format("R = %d Memory usage: %,d bytes %n", numResampledMatrices, (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
        }

        /**
         * Compare against a similar number of Jama.Matrix copies
         */
        System.out.println("Jama.Matrix copies");
        for (int numResampledMatrices = 100; numResampledMatrices <= 100 * 100; numResampledMatrices *= 10) {
            ArrayList<Jama.Matrix> resampledList = new ArrayList<>();

            /**
             * Now, create a resampled data set using that baseline data as the
             * underlying foundation.
             */
            for (int r = 0; r < numResampledMatrices; r++) {
                resampledList.add(baseline.copy());
            }
            System.out.format("R = %d Memory usage: %,d bytes %n", numResampledMatrices, (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
        }

        /**
         * Clearly a massive memory improvement (and speed improvement) over
         * using copies of Jama.Matrix
         */
    }

}
