/*
 * The MIT License
 *
 * Copyright 2018 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class DisjointSet {

    public static void main(String[] args) {
        ArrayList<Integer> data = new ArrayList<>();
        data.add(1);
        data.add(2);
        data.add(3);
        data.add(4);
        // data.add(5);
        // data.add(6);
        // data.add(7);
        // data.add(8);
//        data.add(7);
//        data.add(8);
//        data.add(9);
//        data.add(10);

        // nice, this works
        ArrayList<IntegerPartition> partition = computeAllDisjointSets(data);
        System.out.println(partition);
    }

    public static ArrayList<IntegerPartition> computeAllDisjointSets(ArrayList<Integer> parent) {
        Jama.Matrix binaryCounter = new Jama.Matrix(parent.size(), 1, 0);

        ArrayList<Jama.Matrix> binaryEnumeration = new ArrayList<>();
        count(binaryCounter, 0, binaryEnumeration);

        ArrayList<IntegerPartition> results = new ArrayList<>();
        TreeSet<String> subsetTree = new TreeSet<>();

        for (Jama.Matrix m : binaryEnumeration) {

            ArrayList<Integer> group1 = new ArrayList<>();
            ArrayList<Integer> group2 = new ArrayList<>();
            for (int i = 0; i < m.getRowDimension(); i++) {
                if (m.get(i, 0) == 1) {
                    group1.add(parent.get(i));
                } else {
                    group2.add(parent.get(i));
                }
            }

            String s1 = "{ ";
            for (int i = 0; i < group1.size(); i++) {
                s1 = s1.concat(group1.get(i) + " ");
            }
            s1 = s1.concat("}");
            String s2 = "{ ";
            for (int i = 0; i < group2.size(); i++) {
                s2 = s2.concat(group2.get(i) + " ");
            }
            s2 = s2.concat("}");
            
            boolean debug = false;
            if (debug) {
                System.out.println(s1 + " " + s2);
            }
            IntegerPartition p = new IntegerPartition(group1, group2);

            if (!subsetTree.contains(s1)) {
                if (!subsetTree.contains(s2)) {
                    results.add(p);
                    subsetTree.add(s1);
                    subsetTree.add(s2);
                }
            }
        }
        return results;
    }

    public ArrayList<ArrayList<Integer>> getAllDisjointSets() {
        ArrayList<ArrayList<Integer>> disjoint = new ArrayList<>();

        return disjoint;
    }

    static void count(Matrix binaryCounter, int index, ArrayList<Jama.Matrix> results) {
        // System.out.println("index1: " + index);
        Jama.Matrix ones = new Jama.Matrix(1, binaryCounter.getRowDimension(), 1);

        if (index < binaryCounter.getRowDimension() - 1) {
            count(binaryCounter.copy(), index + 1, results);
            binaryCounter.set(index, 0, 1);
            count(binaryCounter.copy(), index + 1, results);
        } else {
            int maxIndex = (int) Math.floor(binaryCounter.getRowDimension() / 2);
            if ((binaryCounter.getRowDimension() / 2.0) != Math.floor(binaryCounter.getRowDimension() / 2)) {
                maxIndex++;
            }
//            System.out.println(maxIndex);
            // pmUtility.prettyPrintVector(binaryCounter);
            int count = (int) (ones.times(binaryCounter)).get(0, 0);
            boolean valid = true;
            if (count > maxIndex || count < 1) {
                valid = false;
            }
            if (count == maxIndex) {
//                System.out.println("---");
//                pmUtility.prettyPrintVector(binaryCounter);
                boolean repeating = true;
                for (int i = 0; i < 1; i++) {
                    if (binaryCounter.get(i, 0) == 1) {
                        repeating = false;
                    }
                }
                if (repeating) {
                    valid = false;
                }
            }

            if (valid) {
                // have to come up with valid criterion here...
//                System.out.println("*");
                results.add(binaryCounter.copy());
            } else {
//                System.out.println("Rejected");
            }
//            System.out.println("---");

            binaryCounter.set(index, 0, 1);
            // pmUtility.prettyPrintVector(binaryCounter);

            // rule:
            // any one
            // any two, as long as one is one
            // any three, as long as first two
            // any four, as long as first three
            // etc.
            count = (int) (ones.times(binaryCounter)).get(0, 0);
            valid = true;
            if (count > maxIndex || count < 1) {
                valid = false;
            }
            if (count == maxIndex) {
//                System.out.println("---");
//                pmUtility.prettyPrintVector(binaryCounter);
                boolean repeating = true;
                for (int i = 0; i < 1; i++) {
                    if (binaryCounter.get(i, 0) == 1) {
                        repeating = false;
                    }
                }
                if (repeating) {
                    valid = false;
                }
            }

            if (valid) {
//                System.out.println("*");
                results.add(binaryCounter.copy());
            } else {
//                System.out.println("Rejected");
            }
//            System.out.println("---");
        }
    }

}
