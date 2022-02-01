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

import JSci.maths.statistics.ChiSqrDistribution;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeSet;
import optimization.Fmin;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class TreeMoment {

    TreeMoment parent;
    TreeMoment childLeft;
    TreeMoment childRight;
    Boolean isLeftNode;
    Boolean terminal = true;
    Boolean[] discreteVector;

    private SplitRule rule;
    private int depth;

    MomentSpecification momentSpec;
    Jama.Matrix betaEstimateNode;
    Jama.Matrix varianceMatrix;

    final private DataLens lensGrowingTree;
    private DataLens lensHonest;
    int numHonestXObservations;

    int maxDepth = 100;
    double minProportionEachPartition;
    int minCountEachPartition;
    double improvementThreshold;

    boolean verbose;

    boolean debugOptimization = false;
    private double currentNodeObjectiveFunction;
    private ContainerMoment currentNodeMoment;

    private ArrayList<Integer> indexHomogeneousParameters = new ArrayList<>();

    private boolean testParameterHomogeneity;

    public TreeMoment(TreeMoment parent, MomentSpecification spec, DataLens lensGrowingTree,
            Boolean[] discreteVector, boolean verbose, double minProportionEachPartition,
            int minCountEachPartition, double improvementThreshold, boolean isLeft, int maxDepth,
            DataLens lensHonest, boolean testParameterHomogeneity) {
        this.momentSpec = spec;
        this.parent = parent;
        this.lensHonest = lensHonest;
        this.lensGrowingTree = lensGrowingTree;
        this.discreteVector = discreteVector;
        this.verbose = verbose;
        this.minProportionEachPartition = minProportionEachPartition;
        this.minCountEachPartition = minCountEachPartition;
        this.improvementThreshold = improvementThreshold;
        this.isLeftNode = isLeft;
        this.maxDepth = maxDepth;
        this.testParameterHomogeneity = testParameterHomogeneity;

        if (parent == null) {
            depth = 0;
        } else {
            depth = parent.getDepth() + 1;
        }
    }

    /**
     * This will tell us how many treatment effects are estimated in the
     * subtree.
     *
     * @return
     */
    public int countTerminalNodes() {
        if (terminal) {
            return 1;
        } else {
            return childLeft.countTerminalNodes() + childRight.countTerminalNodes();
        }
    }

    /**
     * Return tree's prediction of treatment effect on the basis of x_i and z_i.
     *
     * @param xi Covariates going into the moment function
     * @param zi Covariates that determine the parameter split
     * @return
     */
    public Double getPredictedY(Jama.Matrix xi, Jama.Matrix zi) {
        if (terminal) {
            // return the treatment effect for this x_i, not the predicted y_i
            return momentSpec.getPredictedY(xi, betaEstimateNode);
        } else if (getRule().isLeft(zi)) {
            return childLeft.getPredictedY(xi, zi);
        } else {
            return childRight.getPredictedY(xi, zi);
        }
    }

    public void setDepth(int depth) {
        this.depth = depth;
    }

    public int getDepth() {
        return depth;
    }

    public TreeMoment getChildLeft() {
        return childLeft;
    }

    public TreeMoment getChildRight() {
        return childRight;
    }

    public void setChildLeft(TreeMoment childLeft) {
        this.childLeft = childLeft;
    }

    public void setChildRight(TreeMoment childRight) {
        this.childRight = childRight;
    }

    public void setTerminal(Boolean terminal) {
        this.terminal = terminal;
    }

    //This method builds the tree
    public void determineSplit() {
        if (verbose) {
            echoLn("----------- computing optimal split --------------");
            echoLn(depth + " " + getParentRuleDescriptive(null));
        }

        /**
         * Now compute a baseline SSE for comparing the improvement in fit
         */
        // System.out.println("Computing baseline SSE");
        currentNodeMoment = momentSpec.computeOptimalBeta(lensGrowingTree);
        setNodeEstimatedBeta(currentNodeMoment.getBeta());
        setNodeEstimatedVariance(currentNodeMoment.getVariance());

//        System.out.println("in this sub-node:");
//        if (nodeY.getNumObs() > 10) {
//            pmUtility.prettyPrint(pmUtility.concatMatrix(nodeY, nodeX).getSubsetData(0, 10, 0, nodeX.getColumnDimension()));
//            System.out.println("TreeMoment.java:152 -> min: "+pmUtility.min(nodeX, 1)+" max: "+pmUtility.max(nodeX, 1));
//        }
        /**
         * This is the place to set a priori conditions on growing the tree (max
         * depth, etc.)
         */
        double optimalZ_SSE = 0;
        if (depth < maxDepth) {
            double optimalZ = 0;

            double optimalZ_SSE_Left = 0;
            double optimalZ_SSE_Right = 0;
            int numObsLeft = 0;
            int numObsRight = 0;
            int optimalSplitVariableIndex = 0;
            int optimalDiscreteCollectionIndex = 0;

            /**
             * Go through all the discrete variables and compute all the
             * remaining groups to search over.
             *
             * This is necessary since once you split a discrete variable, you
             * have subgroups that you need to enumerate separately within each
             * node.
             */
            ArrayList<ArrayList<Integer>> discreteCollection = new ArrayList<>();
            ArrayList<Integer> discreteCollectionIndex = new ArrayList<>();
            // for (int k = 0; k < nodeX.getColumnDimension(); k++) {
            for (int k : momentSpec.getVariableIndicesToSearchOver()) { //The variables Z which we can split on
                if (discreteVector[k] == true) {
                    TreeSet<Integer> discreteTreeSet = new TreeSet<>();
                    for (int i = 0; i < lensGrowingTree.getNumObs(); i++) {
                        discreteTreeSet.add((int) lensGrowingTree.getZ(i, k));
                    }
                    ArrayList<Integer> discreteList = new ArrayList<>(discreteTreeSet);
                    discreteCollection.add(discreteList); // discreteCollection is the columns of Z which are discrete
                    // System.out.println("Added list for k = "+k+" of size "+discreteList.size());
                    discreteCollectionIndex.add(k); // discreteCollectionIndex is the indices of the discrete variables
                }
            }

            /**
             * We are going to store p integers here that we will search over.
             * This is going to be random at each split, we don't have to keep
             * track across branches.
             */
            TreeSet<Integer> randomForestIndex = new TreeSet<>();
            boolean useRandomForest = false;
            if (useRandomForest) { // Randomly choosing a subset of variables to search over with a max of 5 variables to search over
                int P = Math.min(5, momentSpec.getVariableIndicesToSearchOver().length);
                for (int i = 0; i < P; i++) {
                    int index = (int) Math.floor(Math.random() * momentSpec.getVariableIndicesToSearchOver().length);
                    while (randomForestIndex.contains(index)) {
                        index = (int) Math.floor(Math.random() * momentSpec.getVariableIndicesToSearchOver().length);
                    }
                    randomForestIndex.add(index);
                }
            } else {
                for (int index : momentSpec.getVariableIndicesToSearchOver()) {
                    randomForestIndex.add(index);
                }
            }

            boolean first = true;

            for (int indexSplitVariable : momentSpec.getVariableIndicesToSearchOver()) {
                if (debugOptimization) {
                    echoLn("indexSplitVariable: " + indexSplitVariable + "; isDiscrete: " + discreteVector[indexSplitVariable] + "; In Tree: " + randomForestIndex.contains(indexSplitVariable));
                }
                if (randomForestIndex.contains(indexSplitVariable)) {
                    /**
                     * Search using fmin and grid search if splitting variable
                     * is continuous
                     */
                    if (discreteVector[indexSplitVariable] == false) {
                        MomentContinuousSplitObj obj = momentSpec.getFminObjective(lensGrowingTree, indexSplitVariable, minProportionEachPartition, minCountEachPartition);
                        double minZ = lensGrowingTree.getMinimumValue(indexSplitVariable);
                        double maxZ = lensGrowingTree.getMaximumValue(indexSplitVariable);
                        // System.out.println("TreeMoment.java:224 -> min x_1: "+minX+" max x_1: "+maxX);
                        double optimalZ_k;
                        double optimalZ_SSE_k;

                        boolean useFmin = true;
                        if (useFmin) {
                            optimalZ_k = Fmin.fmin(minZ, maxZ, obj, 1E-8); // This is choosing a split point such that the summed SSEs of each leaf are minimized
                            optimalZ_SSE_k = obj.f_to_minimize(optimalZ_k); //Now return the summed SSEs of the optimal split point
                        } else {
                            optimalZ_k = Double.POSITIVE_INFINITY;
                            optimalZ_SSE_k = Double.POSITIVE_INFINITY;
                        }

                        if (debugOptimization) {
                            echoLn("\tFmin search on z_" + indexSplitVariable + " found x = " + optimalZ_k + " SSE: " + optimalZ_SSE_k);
                        }

                        boolean testGridSearch = true;
                        double h = 1E-30;
                        if (testGridSearch) {
                            double leftZ = minZ;
                            double rightZ = maxZ;

                            double increment = h + (rightZ - leftZ) / 100.0;

                            if (debugOptimization) {
                                echoLn("\tGrid Search " + momentSpec.getVariableName(indexSplitVariable) + " from " + leftZ + " to " + rightZ);
                            }

                            for (double z = leftZ; z <= rightZ; z += increment) {
                                double f = obj.f_to_minimize(z);
                                if (debugOptimization) {
                                    echoLn("\tGrid search z_" + indexSplitVariable + " (" + momentSpec.getVariableName(indexSplitVariable) + ") from " + optimalZ_SSE_k + " to " + f + " by moving from " + optimalZ_k + " to " + z);
                                }
                                if (f < optimalZ_SSE_k) {
                                    optimalZ_k = z;
                                    optimalZ_SSE_k = f;
                                    // System.out.println("SSE left: "+obj.leftMSE+" SSE right: "+obj.rightMSE);
                                }
                            }
                            if (debugOptimization) {
                                echoLn("\tGrid Search " + momentSpec.getVariableName(indexSplitVariable) + " best " + optimalZ_k + " SSE: " + optimalZ_SSE_k);
                            }

                            /**
                             * Try a second grid search within the last interval
                             * to really improve precision of our estimates
                             */
                            leftZ = optimalZ_k - increment;
                            rightZ = optimalZ_k + increment;
                            increment = h + (rightZ - leftZ) / 100.0;

                            if (debugOptimization) {
                                echoLn("\tSecond Finer Grid Search " + momentSpec.getVariableName(indexSplitVariable) + " from " + leftZ + " to " + rightZ);
                            }

                            for (double z = leftZ; z <= rightZ; z += increment) {
                                double f = obj.f_to_minimize(z);
                                if (debugOptimization) {
                                    echoLn("\tGrid search z_" + indexSplitVariable + " (" + momentSpec.getVariableName(indexSplitVariable) + ") from " + optimalZ_SSE_k + " to " + f + " by moving from " + optimalZ_k + " to " + z);
                                }
                                if (f < optimalZ_SSE_k) {
                                    optimalZ_k = z;
                                    optimalZ_SSE_k = f;
                                    // System.out.println("SSE left: "+obj.leftMSE+" SSE right: "+obj.rightMSE);
                                }
                            }
                            if (debugOptimization) {
                                echoLn("\tGrid Search " + momentSpec.getVariableName(indexSplitVariable) + " best " + optimalZ_k + " SSE: " + optimalZ_SSE_k);
                            }

                        }

                        //If the summed SSE for this variable is smaller than for any other previous variable, or if its the first variable being tested, set it to be the optimal splitting variable
                        if (optimalZ_SSE_k < optimalZ_SSE || first) {
                            obj.f_to_minimize(optimalZ_k);
                            optimalZ = optimalZ_k;
                            optimalZ_SSE = optimalZ_SSE_k;
                            optimalSplitVariableIndex = indexSplitVariable;
                            optimalZ_SSE_Left = obj.getLeftSSE();
                            optimalZ_SSE_Right = obj.getRightSSE();
                            // System.out.println("SSE left: "+optimalZ_SSE_Left+" SSE right: "+optimalZ_SSE_Right);
                            numObsLeft = obj.getEffectiveNumObsLeft();
                            numObsRight = obj.getEffectiveNumObsRight();
                            first = false;
                        }
                    } else {
                        int collectionIndex = discreteCollectionIndex.indexOf(indexSplitVariable);
                        ArrayList<Integer> discreteList = discreteCollection.get(collectionIndex);
                        ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList);
                        // echoLn("Partition size: " + partitions.size() + " discreteList.size(): " + discreteList.size());
                        /**
                         * Need to put in a check here that the
                         * discreteList.size is greater than one element
                         */
                        if (discreteList.size() > 1) {
                            int optimalPartitionIndex = 0;
                            double bestPartitionSSE = 0;
                            double optimalZ_SSE_Right_Partition = 0;
                            double optimalZ_SSE_Left_Partition = 0;
                            int numObsLeft_Partition = 0;
                            int numObsRight_Partition = 0;
                            for (int i = 0; i < partitions.size(); i++) {
                                MomentPartitionObj obj = momentSpec.getMomentPartitionObj(lensGrowingTree, indexSplitVariable, partitions.get(i));

                                double partitionSSE = 0;
                                if (obj.getEffectiveNumObsLeft() < minCountEachPartition || obj.getEffectiveNumObsRight() < minCountEachPartition) {
                                    // echoLn("IS IT IN? : obj.getNumObsLeft(): " + obj.getEffectiveNumObsLeft() + " minCountEachPartition " + minCountEachPartition + " right obs: " + obj.getEffectiveNumObsRight() + " indexSplitVariable " + indexSplitVariable);
                                    if (debugOptimization) {
                                        //   echoLn("\t\tMin K violated: rejecting partition for left obs: " + obj.getNumObsLeft() + " right obs: " + obj.getNumObsRight());
                                    }
                                    partitionSSE = Double.POSITIVE_INFINITY;
                                } else if (((obj.getEffectiveNumObsLeft() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) < minProportionEachPartition || ((obj.getEffectiveNumObsRight() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) < minProportionEachPartition) {
                                    if (debugOptimization) {
                                        //  echoLn("\t\tRejecting partition for proportion; left: " + ((obj.getNumObsLeft() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) + " right: " + ((obj.getNumObsRight() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)));
                                        // System.exit(0);
                                    }
                                    partitionSSE = Double.POSITIVE_INFINITY;
                                } else {
                                    partitionSSE = obj.getSSE();
                                }

                                if (debugOptimization) {
                                    //  echoLn("\t x_" + indexSplitVariable + " Partition: " + i + " " + partitions.get(i) + " SSE: " + partitionSSE);
                                }

                                //For every possible partition, we check which has the lowest SSE
                                if (partitionSSE < bestPartitionSSE || i == 0) {
                                    bestPartitionSSE = partitionSSE;
                                    optimalPartitionIndex = i;
                                    optimalZ_SSE_Right_Partition = obj.getRightSSE();
                                    optimalZ_SSE_Left_Partition = obj.getLeftSSE();
                                    numObsLeft_Partition = obj.getEffectiveNumObsLeft();
                                    numObsRight_Partition = obj.getEffectiveNumObsRight();
                                    if (debugOptimization) {
                                        echoLn("\tPartition: " + i + " " + partitions.get(i) + " SSE: " + partitionSSE + " set as within-variable current best.");
                                    }
                                }
                            }

                            if (bestPartitionSSE < optimalZ_SSE || first) {
                                optimalZ = optimalPartitionIndex; // in this case, this will be the index of the best partitioning
                                // But how do we know which partition the index corresponds to?
                                optimalZ_SSE = bestPartitionSSE;
                                optimalSplitVariableIndex = indexSplitVariable;
                                optimalDiscreteCollectionIndex = collectionIndex;
                                optimalZ_SSE_Left = optimalZ_SSE_Left_Partition;
                                optimalZ_SSE_Right = optimalZ_SSE_Right_Partition;
                                numObsLeft = numObsLeft_Partition;
                                numObsRight = numObsRight_Partition;
                                first = false;
                                if (debugOptimization) {
                                    echoLn("Variable x_" + optimalSplitVariableIndex + " with partition " + partitions.get((int) optimalZ) + " giving SSE of " + optimalZ_SSE + " set as overall current best.");
                                }
                            }
                        }
                    }
                }
            }

            setCurrentNodeObjectiveFunction(currentNodeMoment.getGoodnessOfFit());

            // System.out.println("Baseline SSE is computed as: " + baseline);
            if (verbose) {
                echoLn("Number of observations in node: " + lensGrowingTree.getNumObs());
                echoLn("Improvement from " + getCurrentNodeObjectiveFunction() + " to " + optimalZ_SSE + " (left: " + optimalZ_SSE_Left + " [" + numObsLeft + "] right: " + optimalZ_SSE_Right + " [" + numObsRight + "])");
                echoLn("Improvement absolute (to compare against threshold): " + (getCurrentNodeObjectiveFunction() - optimalZ_SSE));
            }

            // if (baseline - optimalX_SSE < improvementThreshold) {
            if ((getCurrentNodeObjectiveFunction() - optimalZ_SSE) < improvementThreshold || first || getCurrentNodeObjectiveFunction() == 0) {
                setTerminal(true);
                if (verbose) {
                    echoLn(depth + ". Terminating due to lack of improvement in SSE; rules: " + getParentRuleDescriptive(null) + " beta: " + pmUtility.stringPrettyPrintVector(betaEstimateNode));
                }
            } else {
                setTerminal(false);
                if (discreteVector[optimalSplitVariableIndex]) {
                    ArrayList<Integer> discreteList = discreteCollection.get(optimalDiscreteCollectionIndex);
                    ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList); //Recompute all the disjoint sets of the optimal splitting variable, which is discrete

                    MomentPartitionObj obj = momentSpec.getMomentPartitionObj(lensGrowingTree, optimalSplitVariableIndex, partitions.get((int) optimalZ)); //This is where we use optimalX
                    if (verbose) {
                        echoLn(depth + ". Calculated optimal split along discrete variable, partitioning x_" + optimalSplitVariableIndex + " -> " + partitions.get((int) optimalZ) + ", generating SSE of " + obj.getSSE());
                    }
                    setRule(new SplitRule(true, optimalSplitVariableIndex, optimalZ, partitions.get((int) optimalZ), momentSpec));
                    childLeft = new TreeMoment(this, momentSpec, obj.getDataSplit().getLeft(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth, null, testParameterHomogeneity);
                    childRight = new TreeMoment(this, momentSpec, obj.getDataSplit().getRight(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth, null, testParameterHomogeneity);
                } else {
                    MomentContinuousSplitObj obj = momentSpec.getFminObjective(lensGrowingTree, optimalSplitVariableIndex, minProportionEachPartition, minCountEachPartition);
                    if (verbose) {
                        echoLn(depth + ". Calculated optimal split along " + momentSpec.getVariableName(optimalSplitVariableIndex) + " at " + optimalZ + ", generating SSE of " + obj.f_to_minimize(optimalZ));
                    }
                    setRule(new SplitRule(false, optimalSplitVariableIndex, optimalZ, null, momentSpec));

                    DataLens left = SplitContainer.getContinuousDataSplit(lensGrowingTree, optimalZ, optimalSplitVariableIndex).getLeft();
                    DataLens right = SplitContainer.getContinuousDataSplit(lensGrowingTree, optimalZ, optimalSplitVariableIndex).getRight();
                    // System.out.println("Max left: "+pmUtility.max(childLeftX, 1));
                    // System.out.println("Min right: "+pmUtility.min(childRightX, 1));
                    childLeft = new TreeMoment(this, momentSpec, left, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth, null, testParameterHomogeneity);
                    childRight = new TreeMoment(this, momentSpec, right, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth, null, testParameterHomogeneity);
                }
                childLeft.determineSplit(); //Prioritizes splits to the left
                childRight.determineSplit(); //Once we find a terminal node, split to the right
            }
        } else {
            setTerminal(true);
            // System.out.print("before call: ");
            // pmUtility.prettyPrintVector(betaEstimateNode);

            /**
             * Already do all this above
             */
            // ContainerMoment currentNodeMoment = momentSpec.computeOptimalBeta(lensGrowingTree);
            // setNodeEstimatedBeta(currentNodeMoment.getBeta());
            // System.out.print("after call: ");
            // pmUtility.prettyPrintVector(betaEstimateNode);
            // setNodeEstimatedVariance(currentNodeMoment.getVariance());
            if (verbose) {
                echoLn(depth + ". Terminal beta: " + pmUtility.stringPrettyPrintVector(betaEstimateNode));
                // System.out.println(depth + ". Terminal RDD value: " + getRDDEstimate());
            }
        }

    }

//    public String getParentRule(TreeSet<Integer> indexPreviousSplits) {
//        if (terminal && parent == null) {
//            return "{ Stump }";
//        }
//        if (parent == null) {
//            return " }";
//        }
//
//        SplitRule r = parent.getRule();
//        String s = "";
//        if (terminal) {
//            s = "{ ";
//        }
//        if (r.isSplitOnDiscreteVariable()) {
//            if (isLeftNode) {
//                s = s + r.getLeftSplit();
//            } else {
//                s = s + r.getRightSplit();
//            }
//        } else if (isLeftNode) {
//            s = s + "z_" + r.getOptimalSplitVariableIndex() + " < " + r.getSplitPoint();
//        } else {
//            s = s + "z_" + r.getOptimalSplitVariableIndex() + " > " + r.getSplitPoint();
//        }
//        if (indexPreviousSplits == null) {
//            indexPreviousSplits = new TreeSet<>();
//        }
//        if (indexPreviousSplits.contains(r.getOptimalSplitVariableIndex())) { //Don't get this part; this is to avoid reprinting the same rule over and over again for further splits in the tree
//            s = "";
//        } else if (parent.parent != null) {
//            s = s + ", ";
//        }
//
//        indexPreviousSplits.add(r.getOptimalSplitVariableIndex());
//        return s + parent.getParentRule(indexPreviousSplits);
//    }
    public String getParentRuleDescriptive(TreeSet<Integer> indexPreviousSplits) {
        if (terminal && parent == null) {
            return "{ Stump }";
        }
        if (parent == null) {
            return " }";
        }

        SplitRule r = parent.getRule();
        // System.out.println("this node's parent rule: " + r);
        String s = "";
        if (terminal) {
            s = "{ ";
        }
        if (r.isSplitOnDiscreteVariable()) {
            if (isLeftNode) {
                s = s + r.getLeftSplitDescriptive(momentSpec);
            } else {
                s = s + r.getRightSplitDescriptive(momentSpec);
            }
        } else if (isLeftNode) {
            s = s + momentSpec.getVariableName(r.getOptimalSplitVariableIndex()) + " < " + r.getSplitPoint();
        } else {
            s = s + momentSpec.getVariableName(r.getOptimalSplitVariableIndex()) + " > " + r.getSplitPoint();
        }
        if (indexPreviousSplits == null) {
            indexPreviousSplits = new TreeSet<>();
        }
        if (indexPreviousSplits.contains(r.getOptimalSplitVariableIndex())) {
            s = "";
        } else if (parent.parent != null) {
            s = s + ", ";
        }

        indexPreviousSplits.add(r.getOptimalSplitVariableIndex());
        //System.out.println("Current s = " + s + " parent descriptive: " + parent.getParentRuleDescriptive(indexPreviousSplits));
        return s + parent.getParentRuleDescriptive(indexPreviousSplits);
    }

    public void setNodeEstimatedBeta(Jama.Matrix beta) {
        this.betaEstimateNode = beta;
    }

    public Jama.Matrix getNodeEstimatedBeta() {
        return betaEstimateNode;
    }

    public Jama.Matrix getNodeEstimatedBetaVariance() {
        return varianceMatrix;
    }

    public void addAllRulesRecursively(ArrayList<SplitRule> ruleList) {
        if (!terminal) {
            ruleList.add(getRule());
            childLeft.addAllRulesRecursively(ruleList);
            childRight.addAllRulesRecursively(ruleList);
        }
    }

    public void printTree() {
        if (!terminal) {
            childLeft.printTree();
            childRight.printTree();
        } else {
            // System.out.println(getParentRuleDescriptive(null) + " " + pmUtility.stringPrettyPrintVector(betaEstimateNode) + " (" + pmUtility.stringPrettyPrintVector(varianceMatrix) + ")");
            echoLn(getParentRuleDescriptive(null) + " [" + lensHonest.getNumObs() + "] " + momentSpec.formatTreeLeafOutput(betaEstimateNode, varianceMatrix));
        }
    }

    public void distributeHonestObservations(DataLens honest) {
        lensHonest = honest;
        if (!terminal) {
            DataLens[] split = honest.splitOnRule(getRule());
            childLeft.distributeHonestObservations(split[0]);
            childRight.distributeHonestObservations(split[1]);
        }
    }

    public void estimateHonestTree() {
        /**
         * At the top of the tree, start the sorting of honest observations to
         * branches. DataLens should just be able to split on the basis of the
         * rule!
         */
        if (parent == null) {
            distributeHonestObservations(lensHonest);
        }

        if (terminal) {
            if (verbose) {
                echo(getParentRuleDescriptive(null) + " honest set to ");
            }
            if (lensHonest == null) {
                // lensHonest.getNumObs() == 0 || lensHonest ==  null || lensHonest.getXsum(lensHonest.getNumObs()) == 0 || lensHonest.getXsum(lensHonest.getNumObs()) == lensHonest.getNumObs() ) {
                setNodeEstimatedBeta(null);
                setNodeEstimatedVariance(null);
                if (verbose) {
                    echoLn("honestY null");
                }
            } else {
                ContainerMoment c = momentSpec.computeOptimalBeta(lensHonest);
                Jama.Matrix oldBeta = getNodeEstimatedBeta(); //Not sure I understand why this is here... we don't have a beta estimate until we reach a terminal node? A. this is what it was the tree building sample, not the honest tree sample
                setNodeEstimatedBeta(c.getBeta());
                setNodeEstimatedVariance(c.getVariance());

                // for (int i = 0; i < lensHonest.getNumObs(); i++) {
                // echoLn("Did you stop here? 2-6." + " num: " + lensHonest.getNumObs() + " y: " + lensHonest.getY(0) + " x: " + lensHonest.getX(i,0) ); //  + "    " + pmUtility.stringPrettyPrint(c.getBeta()) + "    " + pmUtility.stringPrettyPrint(oldBeta) );
                // }
                if (verbose) {
                    echoLn(pmUtility.stringPrettyPrint(c.getBeta().transpose()) + " [ " + lensHonest.getNumObs() + " ] from " + pmUtility.stringPrettyPrint(oldBeta.transpose()));
                }
            }
        } else {
            childLeft.estimateHonestTree(); //This will climb us down the tree to the terminal nodes
            childRight.estimateHonestTree();
            /**
             * This is where to impose the pruning if the two child nodes have
             * Infinite or null estimates of the treatment effect.
             */
            boolean pruneTree = true;
            if (pruneTree) {
                Jama.Matrix leftBeta = childLeft.getNodeEstimatedBeta();
                Jama.Matrix rightBeta = childRight.getNodeEstimatedBeta();
                if (leftBeta == null || rightBeta == null) {
                    if (leftBeta == null) {
                        if (verbose) {
                            echoLn("Detected child left is null");
                        }
                    }
                    if (rightBeta == null) {
                        if (verbose) {
                            echoLn("Detected child right is null");
                        }
                    }
                    setTerminal(true);
                    if (lensHonest == null) {
                        // lensHonest.getNumObs() == 0 | lensHonest ==  null || lensHonest.getXsum(lensHonest.getNumObs()) == 0 || lensHonest.getXsum(lensHonest.getNumObs()) == lensHonest.getNumObs()     }) {
                        setNodeEstimatedBeta(null);
                        setNodeEstimatedVariance(null);
                        echoLn("null pruned");
                    } else {
                        ContainerMoment c = momentSpec.computeOptimalBeta(lensHonest);
                        setNodeEstimatedBeta(c.getBeta());
                        setNodeEstimatedVariance(c.getVariance());
                        if (verbose) {
                            echoLn("pruned children, replaced with: n = " + lensHonest.getNumObs() + " " + pmUtility.stringPrettyPrintVector(c.getBeta()));
                        }
                    }
                } else {
                    // System.out.println("not pruned: n = " + honestX.getNumObs());
                }
            }
        }

        if (testParameterHomogeneity) {
// System.out.println("DEBUG: Entering testParameterHomogeneity");
            /**
             * Want to think about testing for parameter equality across splits,
             * potentially imposing that homogeneity here (and maybe all nodes
             * below this level?)
             *
             * Ok, I think what we are going to do is do a global imposition via
             * a nested fixed point approach. Outer loop searches over fixed
             * subvector of parameters. Inside loop has the usual tree. We are
             * going to partial out (i.e. subtract it out) the global part, just
             * use the tree for the part where there is (may be) parameter
             * heterogeneity.
             */
            /**
             * For the first split (when the parent is null), test parameter by
             * parameter for homogeneity, then report that
             */
            // System.out.println("Right before Wald test");
            if (parent == null && childLeft != null && childRight != null && 1 == 1) {
                // optimalZ_sse is objective value with heterogeneity
                // test using DM from Newey-McFadden, chi-squared with one degree of freedom
                ChiSqrDistribution chi = new ChiSqrDistribution(1);

                ArrayList<PValue> pList = new ArrayList<>();

                for (int k = 0; k < getNodeEstimatedBeta().getRowDimension(); k++) {
                    DistanceMetricTest wald = new DistanceMetricTest(childLeft.lensHonest.getX(), childRight.lensHonest.getX(), childLeft.lensHonest.getY(), childRight.lensHonest.getY());
                    double dm2 = wald.computeStatistic(k);
                    double pval = 1.0 - chi.cumulative(dm2);
                    System.out.println("p-value: " + pval);
                    pList.add(new PValue(k, pval));
                    if (dm2 < chi.inverse(0.95)) {
                        System.out.println("Absent multiple testing correction, potential parameter homogeneity detected on k = " + k);
                    }

                }
                // now sort the p-values from lowest to highest
                Collections.sort(pList);

                // holm-bonferroni procedure below
                boolean useHolmBonferroni = true;
                if (useHolmBonferroni) {
                    boolean addSuccessiveParameters = false; // need this since i kind of constructed this loop backwards in terms of testing the null hypothesis (which is that there is parameter homogeneity)
                    for (int k = 0; k < pList.size(); k++) {
                        PValue d = pList.get(k);
                        System.out.println(d);
                        double adjustedPValue = 0.05 / (pList.size() - k);
                        System.out.println("p: " + d.getP() + " adjusted P-value: " + adjustedPValue);
                        if (d.getP() > adjustedPValue || addSuccessiveParameters) {
                            System.out.println("Holm-Bonferroni -> Accepting null; Adding parameter index " + d.getK() + " to homogeneity list.");
                            indexHomogeneousParameters.add(d.getK());
                            // should i terminate this for loop here?
                            // i think i am supposed to fail to reject all the other hypotheses if this happens (add them to homogeneity list?)
                            // yes, i sort of wrote this backwards; should be adding parameters to HETEROGENEOUS index
                            addSuccessiveParameters = true;
                        } else {
                            System.out.println("Holm-Bonferroni -> Rejecting null; Retaining parameter index " + d.getK() + " in moment forest.");
                        }
                    }
//                    System.out.println("DEBUG: ending Holm Bonferroni method");
                } else {
                    // easier bonferroni procedure (for checking what's going on here)
                    for (int k = 0; k < pList.size(); k++) {
                        PValue d = pList.get(k);
                        System.out.println(d);
                        double criticalValue = 0.05 / (0.0 + pList.size());
                        System.out.println("p: " + d.getP() + " adjusted critical value: " + criticalValue);
                        if (d.getP() > criticalValue) {
                            System.out.println("Straight Bonferroni -> Accepting null; adding parameter index " + d.getK() + " to homogeneity list.");
                            indexHomogeneousParameters.add(d.getK());
                        } else {
                            System.out.println("Straight Bonferroni -> Rejecting null; retaining parameter index " + d.getK() + " in moment forest.");
                        }
                    }
                }

            }
//            System.out.println("DEBUG: Exiting testParameterHomogeneity");
        }
//        System.out.println("DEBUG: Exiting estimateHonestTree()");
    }

    public void clearEstimationData() {
        /**
         * I am not nulling out the grow lens now because hopefully it is not
         * needed and i am worried about what that does to the master data
         * (probably nothing) but come back to this later if we need to for
         * memory reasons
         */
        if (!terminal) {
            childLeft.clearEstimationData();
            childRight.clearEstimationData();
        }
    }

    public void clearHonestyData() {
        numHonestXObservations = lensHonest.getNumObs();
        if (!terminal) {
            childLeft.clearHonestyData();
            childRight.clearHonestyData();
        }
    }

    public Jama.Matrix getEstimatedBeta(Jama.Matrix zi) {
        if (terminal) {
            // return the treatment effect for this z_i, not the predicted y_i
            return betaEstimateNode;
        } else // use rule to figure out whether to return left or right node's value
        // this will keep going down the rabbit hole until it returns a terminal node's value
        // kind of cool how this works!
        if (getRule().isLeft(zi)) {
            return childLeft.getEstimatedBeta(zi);
        } else {
            return childRight.getEstimatedBeta(zi);
        }
    }

    private void setNodeEstimatedVariance(Jama.Matrix variance) {
        varianceMatrix = variance;
    }

    public Jama.Matrix getVarianceMatrix(Jama.Matrix zi) {
        System.out.println("getVarianceMatrix deprecated");
        System.exit(0);
        if (terminal) {
            // return the treatment effect for this x_i, not the predicted y_i
            return varianceMatrix;
        } else // use rule to figure out whether to return left or right node's value
        // this will keep going down the rabbit hole until it returns a terminal node's value
        // kind of cool how this works!
        if (getRule().isLeft(zi)) {
            return childLeft.getVarianceMatrix(zi);
        } else {
            return childRight.getVarianceMatrix(zi);
        }
    }

    static private void echoLn(String s) {
        // SFIToolkit.displayln(s);
        System.out.println(s);
    }

    static private void echo(String s) {
        // SFIToolkit.displayln(s);
        System.out.print(s);
    }

    /**
     * @return the rule
     */
    public SplitRule getRule() {
        return rule;
    }

    /**
     * @param rule the rule to set
     */
    public void setRule(SplitRule rule) {
        this.rule = rule;
    }

    /**
     * @return the currentNodeObjectiveFunction
     */
    public double getCurrentNodeObjectiveFunction() {
        return currentNodeObjectiveFunction;
    }

    /**
     * @param currentNodeObjectiveFunction the currentNodeObjectiveFunction to
     * set
     */
    public void setCurrentNodeObjectiveFunction(double currentNodeObjectiveFunction) {
        this.currentNodeObjectiveFunction = currentNodeObjectiveFunction;
    }

    /**
     * @return the currentNodeMoment
     */
    public ContainerMoment getCurrentNodeMoment() {
        return currentNodeMoment;
    }

    /**
     * @param currentNodeMoment the currentNodeMoment to set
     */
    public void setCurrentNodeMoment(ContainerMoment currentNodeMoment) {
        this.currentNodeMoment = currentNodeMoment;
    }

    /**
     * @return the indexHomogeneousParameters
     */
    public ArrayList<Integer> getIndexHomogeneousParameters() {
        return indexHomogeneousParameters;
    }

    /**
     * @param testParameterHomogeneity the testParameterHomogeneity to set
     */
    public void setTestParameterHomogeneity(boolean testParameterHomogeneity) {
        this.testParameterHomogeneity = testParameterHomogeneity;
    }
}
