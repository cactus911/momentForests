
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

import core.DataLens;
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
        int numObs = 10;
        Jama.Matrix baselineX = new Jama.Matrix(numObs, 4);
        Jama.Matrix baselineY = new Jama.Matrix(numObs, 1);
        Jama.Matrix balancingVector = new Jama.Matrix(numObs, 1);
        for (int i = 0; i < baselineX.getRowDimension(); i++) {
            for (int j = 0; j < baselineX.getColumnDimension(); j++) {
                baselineX.set(i, j, rng.nextGaussian());
            }
            baselineY.set(i, 0, rng.nextDouble());
            balancingVector.set(i, 0, rng.nextInt(2));
        }

        boolean testSplit = true;
        if(testSplit) {
            testSplit(new DataLens(baselineX, baselineY, balancingVector));
        }
        
        boolean testMemory = false;
        if (testMemory) {
            testMemory(baselineX, baselineY, balancingVector);
        }

        boolean testResampling = true;
        if (testResampling) {
            testResampling(baselineX, baselineY, balancingVector);
        }
    }

    private void testResampling(Jama.Matrix baselineX, Jama.Matrix baselineY, Jama.Matrix balancingVector) {
        /**
         * Show that resampling works
         */
        System.out.println("Original");
        Jama.Matrix concat = pmUtility.concatMatrix(baselineY, baselineX);
        concat = pmUtility.concatMatrix(concat, balancingVector);
        pmUtility.prettyPrint(concat.getMatrix(0, 9, 0, 5));
        System.out.println("Original in DataLens");
        DataLens lens = new DataLens(baselineX, baselineY, balancingVector);
        System.out.println(lens.getSubsetData(0, 9));
        System.out.println("Resampled");
        DataLens resampledLens = lens.getResampledDataLensWithBalance(rng.nextLong());
        System.out.println("Getting subdata");
        DataLens subdata = resampledLens.getSubsetData(0, 9);
        System.out.println("Printing subdata");
        System.out.println(subdata);
        System.exit(0);
    }

    private void testMemory(Jama.Matrix baselineX, Jama.Matrix baselineY, Jama.Matrix balancingVector) {
        /**
         * Compute memory usage using resampled matrices
         */
        System.out.println("ResampledMatrix");
        DataLens lens = new DataLens(baselineX, baselineY, balancingVector);
        
        for (int numResampledMatrices = 100; numResampledMatrices <= 100 * 100; numResampledMatrices *= 10) {
            ArrayList<DataLens> resampledList = new ArrayList<>();

            /**
             * Now, create a resampled data set using that baseline data as the
             * underlying foundation.
             */
            for (int r = 0; r < numResampledMatrices; r++) {
                resampledList.add(lens.getResampledDataLensWithBalance(rng.nextLong()));
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
                resampledList.add(baselineX.copy());
                resampledList.add(baselineY.copy());
                resampledList.add(balancingVector.copy());
            }
            System.out.format("R = %d Memory usage: %,d bytes %n", numResampledMatrices, (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
        }
    }

    private void testSplit(DataLens dataLens) {
        System.out.println("Together");
        System.out.println(dataLens);
        DataLens[] split = dataLens.randomlySplitSample(0.45, rng.nextLong());
        System.out.println("Split 1");
        System.out.println(split[0]);
        System.out.println("Split 2");
        System.out.println(split[1]);
    }

}

