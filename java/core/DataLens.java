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
import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class DataLens {

    final Jama.Matrix originalDataX;
    final Jama.Matrix originalDataY;
    final Jama.Matrix balancingVector;
    int[] dataIndex;

    public DataLens(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix balanceVector) {
        originalDataX = X;
        originalDataY = Y;
        balancingVector = balanceVector;
        dataIndex = new int[originalDataX.getRowDimension()];
        for (int i = 0; i < originalDataX.getRowDimension(); i++) {
            dataIndex[i] = i;
        }
    }

    public DataLens(DataLens d, int[] resampleIndex) {
        originalDataX = d.getOriginalDataX();
        originalDataY = d.getOriginalDataY();
        balancingVector = d.getBalancingVector();
        dataIndex = new int[resampleIndex.length];
        for (int i = 0; i < resampleIndex.length; i++) {
            dataIndex[i] = resampleIndex[i];
        }
    }
    
    public DataLens getDataLensSubset(int[] observations) {
        int[] associatedBackingIndex = new int[observations.length];
        for(int i=0;i<observations.length;i++) {
            associatedBackingIndex[i] = dataIndex[observations[i]];
        }
        return new DataLens(this, associatedBackingIndex);
    }

    private Matrix getOriginalDataX() {
        return originalDataX;
    }

    private Matrix getBalancingVector() {
        return balancingVector;
    }

    private Matrix getOriginalDataY() {
        return originalDataY;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < dataIndex.length; i++) {
            s.append(dataIndex[i] + " ");
            s.append("{ ").append(getY(i)).append(" ");
            for (int j = 0; j < originalDataX.getColumnDimension(); j++) {
                s.append(getX(i, j)).append(" ");
            }
            s.append(balancingVector.get(dataIndex[i], 0)).append(" }\n");
        }
        return s.toString();
    }

    private double getMeanBalancingVector() {
        double ratio = 0;
        // System.out.println("dataIndex.length: "+dataIndex.length);
        for (int i : dataIndex) {
            ratio += balancingVector.get(i, 0);
        }
        ratio /= dataIndex.length;
        // System.out.println(this);
        // System.out.println("balancing ratio: "+ratio); // wait there is something weird here; should always be 0.5
        return ratio;
    }

    public DataLens getResampledDataLensWithBalance(long seed) {
        Random rng = new Random(seed);
        int[] newIndex = new int[dataIndex.length];

        // System.out.println("in getResampledDataLensWithBalance");
        double ratio = getMeanBalancingVector();
        int goalTreatment = (int) Math.round(ratio * dataIndex.length);

        int countTreatment = 0;
        for (int i = 0; i < dataIndex.length; i++) {
            int index = rng.nextInt(dataIndex.length);
            double treatmentIndicator = balancingVector.get(dataIndex[index], 0);
            if (countTreatment < goalTreatment) {
                while (treatmentIndicator == 0) {
                    index = rng.nextInt(dataIndex.length);
                    treatmentIndicator = balancingVector.get(dataIndex[index], 0);
                }
                countTreatment++;
            } else {
                while (treatmentIndicator == 1) {
                    index = rng.nextInt(dataIndex.length);
                    treatmentIndicator = balancingVector.get(dataIndex[index], 0);
                }
            }
            newIndex[i] = index;
        }
        return new DataLens(this, newIndex);
    }

    public DataLens getResampledDataLens(long seed) {
        Random rng = new Random(seed);
        int[] newIndex = new int[dataIndex.length];
        for (int i = 0; i < dataIndex.length; i++) {
            newIndex[i] = rng.nextInt(dataIndex.length);
        }
        return new DataLens(this, newIndex);
    }

    public DataLens[] randomlySplitSampleWithBalance(double proportionFirstSample, long seed) {
        Random rng = new Random(seed);
        DataLens[] splitLens = new DataLens[2]; // the two parts of the split sample
        int sizeFirst = (int) Math.round(proportionFirstSample * dataIndex.length);
        int sizeSecond = dataIndex.length - sizeFirst;

        int[] indicesFirst = new int[sizeFirst];
        int[] indicesSecond = new int[sizeSecond];

        int countTreatmentFirst = 0;
        TreeSet<Integer> treeFirst = new TreeSet<>();

        // System.out.println("in randomlySplitSampleWithBalance");
        int desiredNumTreatmentObsInFirstPart = (int) Math.round(getMeanBalancingVector() * sizeFirst);
        // System.out.println("Aiming for " + desiredNumTreatmentObsInFirstPart + " treatments in first sample.");
        for (int i = 0; i < sizeFirst; i++) {
            int randomInteger = rng.nextInt(getNumObs());
            double treatmentStatus = balancingVector.get(dataIndex[randomInteger], 0);
            if (countTreatmentFirst < desiredNumTreatmentObsInFirstPart) {
                while (treatmentStatus == 0 || treeFirst.contains(randomInteger)) {
                    randomInteger = rng.nextInt(getNumObs());
                    treatmentStatus = balancingVector.get(dataIndex[randomInteger], 0);
                }
                treeFirst.add(randomInteger);
                countTreatmentFirst++;
            } else {
                while (treatmentStatus == 1 || treeFirst.contains(randomInteger)) {
                    randomInteger = rng.nextInt(getNumObs());
                    treatmentStatus = balancingVector.get(dataIndex[randomInteger], 0);
                }
                treeFirst.add(randomInteger);
            }
        }
        ArrayList<Integer> firstList = new ArrayList(treeFirst);
        // System.out.println(firstList.size() + " out of desired " + sizeFirst);
        for (int i = 0; i < sizeFirst; i++) {
            indicesFirst[i] = dataIndex[firstList.get(i)];
        }

        int secondCounter = 0;
        for (int i = 0; i < getNumObs(); i++) {
            if (!treeFirst.contains(i)) {
                indicesSecond[secondCounter] = dataIndex[i];
                secondCounter++;
            }
        }

        splitLens[0] = new DataLens(this, indicesFirst);
        splitLens[1] = new DataLens(this, indicesSecond);

        return splitLens;
    }

    /**
     * Obtain the number of observations in the DataLens.
     *
     * @return Number of observations.
     */
    public int getNumObs() {
        return dataIndex.length;
    }

    public int getColumnDimensionX() {
        // System.out.println("Getting column dimension");
        return originalDataX.getColumnDimension();
    }

    public DataLens getSubsetData(int rowStart, int rowEnd) {
        // i am going to just restrict the index to between rowstart and rowend and then return that lens
        int[] subIndex = new int[rowEnd - rowStart + 1];
        int counter = 0;
        for (int i = rowStart; i <= rowEnd; i++) {
            subIndex[counter] = dataIndex[i];
            counter++;
        }
        return new DataLens(this, subIndex);
    }

    public double getY(int i) {
        return originalDataY.get(dataIndex[i], 0);
    }

    public double getX(int i, int j) {
        return originalDataX.get(dataIndex[i], j);
    }

    public Jama.Matrix getRowX(int row) {
        return originalDataX.getMatrix(dataIndex[row], dataIndex[row], 0, originalDataX.getColumnDimension() - 1);
    }

    double getMinimumValue(int indexSplitVariable) {
        double minimumValue = getX(0, indexSplitVariable);
        for (int i = 1; i < getNumObs(); i++) {
            double v = getX(i, indexSplitVariable);
            if (v < minimumValue) {
                minimumValue = v;
            }
        }
        return minimumValue;
    }

    double getMaximumValue(int indexSplitVariable) {
        double maximumValue = getX(0, indexSplitVariable);
        for (int i = 1; i < getNumObs(); i++) {
            double v = getX(i, indexSplitVariable);
            if (v > maximumValue) {
                maximumValue = v;
            }
        }
        return maximumValue;
    }

    Matrix getRowAsJamaMatrix(int i) {
        Jama.Matrix rowX = new Jama.Matrix(1, originalDataX.getColumnDimension());
        for (int j = 0; j < originalDataX.getColumnDimension(); j++) {
            rowX.set(0, j, getX(i, j));
        }
        return rowX;
    }

    DataLens[] splitOnRule(SplitRule rule) {
        DataLens[] split = new DataLens[2];
        ArrayList<Integer> leftList = new ArrayList<>();
        ArrayList<Integer> rightList = new ArrayList<>();
        for (int i = 0; i < getNumObs(); i++) {
            Jama.Matrix xi = getRowAsJamaMatrix(i);
            if (rule.isLeft(xi)) {
                leftList.add(dataIndex[i]);
            } else {
                rightList.add(dataIndex[i]);
            }
        }
        int[] leftIndices = new int[leftList.size()];
        int[] rightIndices = new int[rightList.size()];
        for (int i = 0; i < leftList.size(); i++) {
            leftIndices[i] = leftList.get(i);
        }
        for (int i = 0; i < rightList.size(); i++) {
            rightIndices[i] = rightList.get(i);
        }
        split[0] = new DataLens(this, leftIndices);
        split[1] = new DataLens(this, rightIndices);
        return split;
    }

}
