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
package core;

import Jama.Matrix;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ResamplingLens extends Jama.Matrix {

    final Jama.Matrix originalData;
    long seed;
    Random rng = new Random();
    int[] resampleIndex;

    public ResamplingLens(Jama.Matrix m, long seed) {
        super(0,0);
        originalData = m; // pass by reference
        rng.setSeed(seed);
        initializeIndex();
    }
    
    public ResamplingLens(Jama.Matrix m, long seed, Jama.Matrix balancing) {
        super(0,0);
        originalData = m; // pass by reference
        rng.setSeed(seed);
        initializeIndex(balancing);
    }
    
    /**
     * Use this constructor to create a submatrix. Simply pass the int[] index to the
     * constructor to get a subset of the original data.
     * 
     * @param m Original data Jama.Matrix.
     * @param seed Random number generator seed.
     * @param resampleIndex An index that controls re/subsampling
     */
    public ResamplingLens(Jama.Matrix m, long seed, int[] resampleIndex) {
        super(0,0);
        originalData = m;
        rng.setSeed(seed); // don't even think this is necessary for the subsampling part
        this.resampleIndex = resampleIndex;
    }

    private void initializeIndex(Jama.Matrix balancingVector) {
        resampleIndex = new int[originalData.getRowDimension()];
        double ratio = pmUtility.mean(balancingVector, 0);
        int goalTreatment = (int) Math.round(ratio * balancingVector.getRowDimension());

        int countTreatment = 0;
        for (int i = 0; i < originalData.getRowDimension(); i++) {
            int index = rng.nextInt(originalData.getRowDimension());
            double treatmentIndicator = balancingVector.get(index, 0);
            if (countTreatment < goalTreatment) {
                while (treatmentIndicator == 0) {
                    index = rng.nextInt(originalData.getRowDimension());
                    treatmentIndicator = balancingVector.get(index, 0);
                }
                countTreatment++;
            } else {
                while (treatmentIndicator == 1) {
                    index = rng.nextInt(originalData.getRowDimension());
                    treatmentIndicator = balancingVector.get(index, 0);
                }
            }
            resampleIndex[i] = index;
        }
    }
    
    private void initializeIndex() {
        resampleIndex = new int[originalData.getRowDimension()];
        for (int i = 0; i < originalData.getRowDimension(); i++) {
            resampleIndex[i] = rng.nextInt(originalData.getRowDimension());
        }
    }

    @Override
    public double get(int row, int column) {
        return originalData.get(resampleIndex[row], column);
    }

    /**
     * We will control submatrices through here. 
     * @return Number of observations.
     */
    @Override
    public int getRowDimension() {
        return resampleIndex.length;
    }

    @Override
    public int getColumnDimension() {
        // System.out.println("Getting column dimension");
        return originalData.getColumnDimension();
    }

    @Override
    public double[][] getArray() {
        System.out.println("getArray");
        System.exit(0);
        return super.getArray();
    }

    @Override
    public double[][] getArrayCopy() {
        System.out.println("getArrayCopy");
        System.exit(0);
        return super.getArrayCopy(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getColumnPackedCopy() {
        System.out.println("getColumnPackedCopy");
        System.exit(0);
        return super.getColumnPackedCopy(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix getMatrix(int[] r, int[] c) {
        System.out.println("getMatrix2");
        System.exit(0);
        return super.getMatrix(r, c); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix getMatrix(int i0, int i1, int[] c) {
        System.out.println("getMatrix3");
        System.exit(0);
        return super.getMatrix(i0, i1, c); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix getMatrix(int[] r, int j0, int j1) {
        System.out.println("getMatrix4");
        System.exit(0);
        return super.getMatrix(r, j0, j1); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getRowPackedCopy() {
        System.out.println("getRowPackedCopy");
        System.exit(0);
        return super.getRowPackedCopy(); //To change body of generated methods, choose Tools | Templates.
    }
    

    @Override
    public Jama.Matrix getMatrix(int rowStart, int rowEnd, int columnStart, int columnEnd) {
        Jama.Matrix submatrix = new Jama.Matrix(rowEnd - rowStart + 1, columnEnd - columnStart + 1);
        int counter = 0;
        for (int i = rowStart; i <= rowEnd; i++) {
            int originalIndex = resampleIndex[i];
            for (int j = columnStart; j <= columnEnd; j++) {
                submatrix.set(counter, j, originalData.get(originalIndex, j));
            }
            counter++;
        }
        return submatrix;
    }

}
