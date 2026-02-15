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
import java.util.Collections;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.Arrays;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class DataLens {

    final Jama.Matrix originalDataX;
    final Jama.Matrix originalDataY;
    final Jama.Matrix originalDataZ;
    final Jama.Matrix balancingVector;
    final int[] strataColumnIndex;
    int[] dataIndex;

    //DataLens for original data
    public DataLens(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z, Jama.Matrix balanceVector) {
        originalDataX = X;
        originalDataY = Y;
        originalDataZ = Z;
        strataColumnIndex = null;
        balancingVector = balanceVector; // Is the balance vector specific to RCT and natural experiment settings?
        dataIndex = new int[originalDataX.getRowDimension()];
        for (int i = 0; i < originalDataX.getRowDimension(); i++) {
            dataIndex[i] = i;
        }
    }

    public DataLens(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z, Jama.Matrix balanceVector, int[] strataColIndex) {
        originalDataX = X;
        originalDataY = Y;
        originalDataZ = Z;
        balancingVector = balanceVector;
        strataColumnIndex = strataColIndex;
        dataIndex = new int[originalDataX.getRowDimension()];
        for (int i = 0; i < originalDataX.getRowDimension(); i++) {
            dataIndex[i] = i;
        }
    }

//    public DataLens(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z) {
//        originalDataX = X;
//        originalDataY = Y;
//        originalDataZ = Z;
//        balancingVector = null;
//        dataIndex = new int[originalDataX.getRowDimension()];
//        for (int i = 0; i < originalDataX.getRowDimension(); i++) {
//            dataIndex[i] = i;
//        }
//    }
    //Datalens for the resampled tree data or for a leaf once split
    //Does each leaf carry around the original data then an index of the observations currently in that leaf?
    public DataLens(DataLens d, int[] resampleIndex) {
        originalDataX = d.getOriginalDataX();
        originalDataY = d.getOriginalDataY();
        originalDataZ = d.getOriginalDataZ();
        balancingVector = d.getBalancingVector();
        strataColumnIndex = d.getStrataColumnIndex();
        dataIndex = new int[resampleIndex.length];
        for (int i = 0; i < resampleIndex.length; i++) {
            dataIndex[i] = resampleIndex[i];
        }
    }

    //Gets a subset of the current datalens indexed by "observations"
    public DataLens getDataLensSubset(int[] observations) {
        int[] associatedBackingIndex = new int[observations.length];
        for (int i = 0; i < observations.length; i++) {
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

    private Matrix getOriginalDataZ() {
        return originalDataZ;
    }

    public int[] getStrataColumnIndex() {
        return strataColumnIndex;
    }

    //Prints the original data Y X B, where B is the balancing vector
    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < dataIndex.length; i++) {
            s.append(dataIndex[i]);
            s.append(" ");
            s.append("{ ").append(getY(i)).append(" ");
            for (int j = 0; j < originalDataX.getColumnDimension(); j++) {
                s.append(getX(i, j)).append(" ");
            }
            for (int j = 0; j < originalDataZ.getColumnDimension(); j++) {
                s.append(getZ(i, j)).append(" ");
            }
            if (balancingVector != null) {
                s.append(balancingVector.get(dataIndex[i], 0)).append(" }\n");
            }
        }
        return s.toString();
    }

    //Returns the fraction of treated obs
    private double getMeanBalancingVector() {
        if (balancingVector == null) {
            throw new IllegalStateException("Querying model for average balancing vector when it is null.");
        }
        double ratio = 0;
        for (int i : dataIndex) {
            ratio += balancingVector.get(i, 0);
        }
        ratio /= dataIndex.length;
        return ratio;
    }

    //Performs the resampling for tree data, keeping the proportion of treated and untreated obs the same as in original data and returns the resulting datalense
    public DataLens getResampledDataLensWithBalance(long seed) {
        if (balancingVector == null) {
            throw new IllegalStateException("Trying to resample with average balancing vector when it is null.");
        }
        Random rng = new Random(seed);
        // int[] newIndex = new int[originalDataX.getRowDimension()];
        int[] newIndex = new int[dataIndex.length];

        double ratio = getMeanBalancingVector();
        int goalTreatment = (int) Math.round(ratio * dataIndex.length);

        // Check that the data contains both treated and untreated observations
        boolean hasTreated = false;
        boolean hasUntreated = false;
        for (int idx : dataIndex) {
            if (balancingVector.get(idx, 0) == 1) {
                hasTreated = true;
            } else {
                hasUntreated = true;
            }
            if (hasTreated && hasUntreated) {
                break;
            }
        }
        if (goalTreatment > 0 && !hasTreated) {
            throw new IllegalStateException("Cannot resample with balance: no treated observations in data but goalTreatment = " + goalTreatment);
        }
        if (goalTreatment < dataIndex.length && !hasUntreated) {
            throw new IllegalStateException("Cannot resample with balance: no untreated observations in data but need " + (dataIndex.length - goalTreatment) + " untreated");
        }

        int countTreatment = 0;
        for (int i = 0; i < dataIndex.length; i++) {
            int index = rng.nextInt(dataIndex.length);
            double treatmentIndicator = balancingVector.get(dataIndex[index], 0);
            if (countTreatment < goalTreatment) {
                while (treatmentIndicator == 0) { // Keep resampling until getting a treated obs
                    index = rng.nextInt(dataIndex.length);
                    treatmentIndicator = balancingVector.get(dataIndex[index], 0);
                }
                countTreatment++;
            } else {
                while (treatmentIndicator == 1) { // Now keep resampling until getting an untreated obs
                    index = rng.nextInt(dataIndex.length);
                    treatmentIndicator = balancingVector.get(dataIndex[index], 0);
                }
            }
            newIndex[i] = index;
        }
        return new DataLens(this, newIndex);
    }

    //Performs the resampling for tree data without balancing treated and untreated; or if there is no treatment variable and returns the resulting datalens
    public DataLens getResampledDataLens(long seed) {
        Random rng = new Random(seed);
        int[] newIndex = new int[dataIndex.length];
        for (int i = 0; i < dataIndex.length; i++) {
            newIndex[i] = dataIndex[rng.nextInt(dataIndex.length)];
        }
        return new DataLens(this, newIndex);
    }

    public DataLens getSubsampledDataLens(long seed, double d) {
        Random rng = new Random(seed);
        int b = (int) Math.round(Math.pow(dataIndex.length, d));
        b = Math.min(b, dataIndex.length);
        int[] newIndex = new int[b];
        TreeSet<Integer> drawnTree = new TreeSet<>();
        for (int i = 0; i < b; i++) {
            int candidate = rng.nextInt(dataIndex.length);
            while (drawnTree.contains(candidate)) {
                candidate = rng.nextInt(dataIndex.length);
            }
            newIndex[i] = dataIndex[candidate];
            drawnTree.add(candidate);
        }
        return new DataLens(this, newIndex);
    }

    //Randomly splits the data into the growing and estimating samples and return them
    //Input is the randomly resampled data
    //Output is a datalens vector containing the datalens for each sample
    public DataLens[] randomlySplitSampleWithBalance(double proportionFirstSample, long seed) {
        Random rng = new Random(seed);
        DataLens[] splitLens = new DataLens[2]; // the two parts of the split sample
        int sizeFirst = (int) Math.round(proportionFirstSample * dataIndex.length); //Number of obs going to growing sample
        int sizeSecond = dataIndex.length - sizeFirst;

        int[] indicesFirst = new int[sizeFirst];
        int[] indicesSecond = new int[sizeSecond];

        int countTreatmentFirst = 0;
        TreeSet<Integer> treeFirst = new TreeSet<>(); //Why do we use this TreeSet?

        int desiredNumTreatmentObsInFirstPart = (int) Math.round(getMeanBalancingVector() * sizeFirst);

        // Check that the data contains enough treated and untreated observations
        int totalTreated = 0;
        int totalUntreated = 0;
        for (int idx : dataIndex) {
            if (balancingVector.get(idx, 0) == 1) {
                totalTreated++;
            } else {
                totalUntreated++;
            }
        }
        if (desiredNumTreatmentObsInFirstPart > 0 && totalTreated == 0) {
            throw new IllegalStateException("Cannot split sample with balance: no treated observations in data but need " + desiredNumTreatmentObsInFirstPart + " treated in first sample");
        }
        int desiredUntreatedFirst = sizeFirst - desiredNumTreatmentObsInFirstPart;
        if (desiredUntreatedFirst > 0 && totalUntreated == 0) {
            throw new IllegalStateException("Cannot split sample with balance: no untreated observations in data but need " + desiredUntreatedFirst + " untreated in first sample");
        }

        // System.out.println("Aiming for " + desiredNumTreatmentObsInFirstPart + " treatments in first sample.");
        for (int i = 0; i < sizeFirst; i++) {
            int guess = rng.nextInt(getNumObs());
            double treatmentStatus = balancingVector.get(dataIndex[guess], 0);
            if (countTreatmentFirst < desiredNumTreatmentObsInFirstPart) {
                while (treatmentStatus == 0 || treeFirst.contains(guess)) {
                    guess = rng.nextInt(getNumObs());
                    treatmentStatus = balancingVector.get(dataIndex[guess], 0);
                }
                treeFirst.add(guess);
                countTreatmentFirst++;
            } else {
                while (treatmentStatus == 1 || treeFirst.contains(guess)) {
                    guess = rng.nextInt(getNumObs());
                    treatmentStatus = balancingVector.get(dataIndex[guess], 0);
                }
                treeFirst.add(guess);
            }
        }
        ArrayList<Integer> firstList = new ArrayList<>(treeFirst);
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

        splitLens[0] = new DataLens(this, indicesFirst); //Creates resampled datalens for growing tree
        splitLens[1] = new DataLens(this, indicesSecond); //Creates resampled datalens for estimating tree

        return splitLens;
    }

    /**
     *
     * Randomly split the sample into two parts depending on parameters below.
     * This method does not balance based on balanceVector.
     *
     * @param proportionFirstSample What proportion of the sample to split into
     * the first chunk.
     * @param seed A random number seed.
     * @return
     */
    public DataLens[] randomlySplitSample(double proportionFirstSample, long seed) {
        Random rng = new Random(seed);
        DataLens[] splitLens = new DataLens[2]; // the two parts of the split sample
        int sizeFirst = (int) Math.round(proportionFirstSample * dataIndex.length); //Number of obs going to first sample
        int sizeSecond = dataIndex.length - sizeFirst;

        int[] indicesFirst = new int[sizeFirst];
        int[] indicesSecond = new int[sizeSecond];

        TreeSet<Integer> treeFirst = new TreeSet<>(); // this keeps track of which indices are in the first sample

        for (int i = 0; i < sizeFirst; i++) {
            int guess = rng.nextInt(getNumObs());
            while (treeFirst.contains(guess)) {
                guess = rng.nextInt(getNumObs());
            }
            treeFirst.add(guess);
        }
        ArrayList<Integer> firstList = new ArrayList<>(treeFirst);
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

        splitLens[0] = new DataLens(this, indicesFirst); //Creates resampled datalens for growing tree
        splitLens[1] = new DataLens(this, indicesSecond); //Creates resampled datalens for estimating tree

        return splitLens;
    }

    /**
     *
     * Randomly split the sample into two parts depending on parameters below.
     * This method performs stratified random sampling based on column indices
     *
     * @param proportionFirstSample What proportion of the sample to split into
     * the first chunk.
     * @param seed A random number seed.
     * @return
     */
    public DataLens[] randomlySplitSampleByStrata(double proportionFirstSample, long seed) {
        if (strataColumnIndex == null || strataColumnIndex.length == 0) {
            throw new IllegalStateException("Strata column indices not set in DataLens.");
        }

        Random rng = new Random(seed);
        ArrayList<Integer> indicesFirstList = new ArrayList<>();
        ArrayList<Integer> indicesSecondList = new ArrayList<>();

        // Group observations by combined strata key
        Map<String, List<Integer>> strataGroups = new HashMap<>();
        for (int obsIndex : dataIndex) {
            String key = buildStrataKey(obsIndex, strataColumnIndex);
            List<Integer> group = strataGroups.get(key);
            if (group == null) {
                group = new ArrayList<>();
                strataGroups.put(key, group);
            }
            group.add(obsIndex);
        }

        // Split within each group
        for (List<Integer> group : strataGroups.values()) {
            Collections.shuffle(group, rng);
            int cutoff = (int) Math.round(group.size() * proportionFirstSample);
            indicesFirstList.addAll(group.subList(0, cutoff));
            indicesSecondList.addAll(group.subList(cutoff, group.size()));
        }

        int[] indicesFirst = indicesFirstList.stream().mapToInt(Integer::intValue).toArray();
        int[] indicesSecond = indicesSecondList.stream().mapToInt(Integer::intValue).toArray();

        DataLens[] splitLens = new DataLens[2];
        splitLens[0] = new DataLens(this, indicesFirst);
        splitLens[1] = new DataLens(this, indicesSecond);
        return splitLens;
    }

    private String buildStrataKey(int obsIndex, int[] strataColumnIndices) {
        StringBuilder sb = new StringBuilder();
        for (int col : strataColumnIndices) {
            sb.append(originalDataZ.get(obsIndex, col)).append("_");
        }
        return sb.toString();
    }

    /**
     * Obtain the number of observations in the DataLens.
     *
     * @return Number of observations.
     */
    public int getNumObs() {
        return dataIndex.length;
    }

    public Jama.Matrix getX() {
        // generate a new matrix using the dataIndex
        Jama.Matrix tempX = new Jama.Matrix(getNumObs(), getColumnDimensionX());
        for (int i = 0; i < getNumObs(); i++) {
            for (int j = 0; j < getColumnDimensionX(); j++) {
                tempX.set(i, j, originalDataX.get(dataIndex[i], j));
            }
        }
        return tempX;
    }

    public Jama.Matrix getZ() {
        // generate a new matrix using the dataIndex
        Jama.Matrix tempZ = new Jama.Matrix(getNumObs(), getColumnDimensionZ());
        for (int i = 0; i < getNumObs(); i++) {
            for (int j = 0; j < getColumnDimensionZ(); j++) {
                tempZ.set(i, j, originalDataZ.get(dataIndex[i], j));
            }
        }
        return tempZ;
    }

    public double getZ(int i, int j) {
        return originalDataZ.get(dataIndex[i], j);
    }

    public Jama.Matrix getY() {
        // generate a new matrix using the dataIndex
        Jama.Matrix tempY = new Jama.Matrix(getNumObs(), 1);
        for (int i = 0; i < getNumObs(); i++) {
            tempY.set(i, 0, originalDataY.get(dataIndex[i], 0));
        }
        return tempY;
    }

    public int getColumnDimensionX() {
        // System.out.println("Getting column dimension");
        return originalDataX.getColumnDimension();
    }

    public int getColumnDimensionZ() {
        // System.out.println("Getting column dimension");
        return originalDataZ.getColumnDimension();
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

    public Jama.Matrix getRowZ(int row) {
        return originalDataZ.getMatrix(dataIndex[row], dataIndex[row], 0, originalDataZ.getColumnDimension() - 1);
    }

    //Returns the minimum value of the split variable
    double getMinimumValue(int indexSplitVariable) {
        double minimumValue = getZ(0, indexSplitVariable);
        for (int i = 1; i < getNumObs(); i++) {
            double v = getZ(i, indexSplitVariable);
            if (v < minimumValue) {
                minimumValue = v;
            }
        }
        return minimumValue;
    }

    double getMaximumValue(int indexSplitVariable) {
        double maximumValue = getZ(0, indexSplitVariable);
        for (int i = 1; i < getNumObs(); i++) {
            double v = getZ(i, indexSplitVariable);
            if (v > maximumValue) {
                maximumValue = v;
            }
        }
        return maximumValue;
    }

    Matrix getRowXAsJamaMatrix(int i) {
        Jama.Matrix rowX = new Jama.Matrix(1, originalDataX.getColumnDimension());
        for (int j = 0; j < originalDataX.getColumnDimension(); j++) {
            rowX.set(0, j, getX(i, j));
        }
        return rowX;
    }

    Matrix getRowZAsJamaMatrix(int i) {
        Jama.Matrix rowZ = new Jama.Matrix(1, originalDataZ.getColumnDimension());
        for (int j = 0; j < originalDataZ.getColumnDimension(); j++) {
            rowZ.set(0, j, getZ(i, j));
        }
        return rowZ;
    }

    //Performs the splitting of the data based on the optimal splitting rule and returns the datalens for each leaf
    DataLens[] splitOnRule(SplitRule rule) {
        DataLens[] split = new DataLens[2];
        ArrayList<Integer> leftList = new ArrayList<>();
        ArrayList<Integer> rightList = new ArrayList<>();
        for (int i = 0; i < getNumObs(); i++) {
            Jama.Matrix zi = getRowZAsJamaMatrix(i);
            if (rule.isLeft(zi)) {
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
