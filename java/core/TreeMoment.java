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

import java.util.ArrayList;
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

    // private final static double improvementThreshold = 0.01;
    int maxDepth = 100;
    double minProportionEachPartition;
    int minCountEachPartition;
    double improvementThreshold;

    boolean verbose;

    boolean debugOptimization = false;

    public TreeMoment(TreeMoment parent, MomentSpecification spec, DataLens lensGrowingTree,
            Boolean[] discreteVector, boolean verbose, double minProportionEachPartition,
            int minCountEachPartition, double improvementThreshold, boolean isLeft, int maxDepth,
            DataLens lensHonest) {
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
         * Now compute a baseline MSE for comparing the improvement in fit
         */
        // System.out.println("Computing baseline SSE");
        ContainerMoment currentNodeMoment = momentSpec.computeOptimalBeta(lensGrowingTree);
        setNodeEstimatedBeta(currentNodeMoment.getBeta());
        setNodeEstimatedVariance(currentNodeMoment.getVariance());

        double baseline = currentNodeMoment.getMSE();

//        System.out.println("in this sub-node:");
//        if (nodeY.getNumObs() > 10) {
//            pmUtility.prettyPrint(pmUtility.concatMatrix(nodeY, nodeX).getSubsetData(0, 10, 0, nodeX.getColumnDimension()));
//            System.out.println("TreeMoment.java:152 -> min: "+pmUtility.min(nodeX, 1)+" max: "+pmUtility.max(nodeX, 1));
//        }
        /**
         * This is the place to set a priori conditions on growing the tree (max
         * depth, etc.)
         */
        if (depth < maxDepth) {
            double optimalZ = 0;
            double optimalZ_MSE = 0;
            double optimalZ_MSE_Left = 0;
            double optimalZ_MSE_Right = 0;
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
                        double optimalZ_MSE_k;

                        boolean useFmin = true;
                        if (useFmin) {
                            optimalZ_k = Fmin.fmin(minZ, maxZ, obj, 1E-100); // This is choosing a split point such that the summed MSEs of each leaf are minimized
                            optimalZ_MSE_k = obj.f_to_minimize(optimalZ_k); //Now return the summed MSEs of the optimal split point
                        } else {
                            optimalZ_k = Double.POSITIVE_INFINITY;
                            optimalZ_MSE_k = Double.POSITIVE_INFINITY;
                        }

                        if (debugOptimization) {
                            echoLn("\tFmin search on z_" + indexSplitVariable + " found x = " + optimalZ_k + " mse: " + optimalZ_MSE_k);
                        }

                        boolean testGridSearch = true;
                        double h = 1E-30;
                        if (testGridSearch) {

                            double increment = h + (maxZ - minZ) / 100.0;
                            
                            if (debugOptimization) {
                                echoLn("\tGrid Search " + momentSpec.getVariableName(indexSplitVariable) + " from " + minZ + " to " + maxZ);
                            }

                            for (double z = minZ; z <= maxZ; z += increment) {
                                double f = obj.f_to_minimize(z);
                                if (debugOptimization) {
                                    echoLn("\tGrid search x_" + indexSplitVariable + " (" + momentSpec.getVariableName(indexSplitVariable) + ") from " + optimalZ_MSE_k + " to " + f + " by moving from " + optimalZ_k + " to " + z);
                                }
                                if (f < optimalZ_MSE_k) {
                                    optimalZ_k = z;
                                    optimalZ_MSE_k = f;
                                }
                            }
                            if (debugOptimization) {
                                echoLn("\tGrid Search " + momentSpec.getVariableName(indexSplitVariable) + " best " + optimalZ_k + " MSE: " + optimalZ_MSE_k);
                            }

                        }

                        //If the summed MSE for this variable is smaller than for any other previous variable, or if its the first variable being tested, set it to be the optimal splitting variable
                        if (optimalZ_MSE_k < optimalZ_MSE || first) {
                            optimalZ = optimalZ_k;
                            optimalZ_MSE = optimalZ_MSE_k;
                            optimalSplitVariableIndex = indexSplitVariable;
                            optimalZ_MSE_Left = obj.getLeftMSE();
                            optimalZ_MSE_Right = obj.getRightMSE();
                            numObsLeft = obj.getEffectiveNumObsLeft();
                            numObsRight = obj.getEffectiveNumObsRight();
                            first = false;
                        }
                    } else {
                        int collectionIndex = discreteCollectionIndex.indexOf(indexSplitVariable);
                        ArrayList<Integer> discreteList = discreteCollection.get(collectionIndex);
                        ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList);
                        echoLn("Partition size: " + partitions.size() + " discreteList.size(): " + discreteList.size());
                        /**
                         * Need to put in a check here that the
                         * discreteList.size is greater than one element
                         */
                        if (discreteList.size() > 1) {
                            int optimalPartitionIndex = 0;
                            double bestPartitionMSE = 0;
                            double optimalZ_MSE_Right_Partition = 0;
                            double optimalZ_MSE_Left_Partition = 0;
                            int numObsLeft_Partition = 0;
                            int numObsRight_Partition = 0;
                            for (int i = 0; i < partitions.size(); i++) {
                                MomentPartitionObj obj = momentSpec.getMomentPartitionObj(lensGrowingTree, indexSplitVariable, partitions.get(i));

                                double partitionMSE = 0;
                                if (obj.getEffectiveNumObsLeft() < minCountEachPartition || obj.getEffectiveNumObsRight() < minCountEachPartition) {
                                    echoLn("IS IT IN? : obj.getNumObsLeft(): " + obj.getEffectiveNumObsLeft() + " minCountEachPartition " + minCountEachPartition + " right obs: " + obj.getEffectiveNumObsRight() + " indexSplitVariable " + indexSplitVariable);
                                    if (debugOptimization) {
                                        //   echoLn("\t\tMin K violated: rejecting partition for left obs: " + obj.getNumObsLeft() + " right obs: " + obj.getNumObsRight());
                                    }
                                    partitionMSE = Double.POSITIVE_INFINITY;
                                } else if (((obj.getEffectiveNumObsLeft() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) < minProportionEachPartition || ((obj.getEffectiveNumObsRight() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) < minProportionEachPartition) {
                                    if (debugOptimization) {
                                        //  echoLn("\t\tRejecting partition for proportion; left: " + ((obj.getNumObsLeft() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)) + " right: " + ((obj.getNumObsRight() + 0.0) / (lensGrowingTree.getNumObs() + 0.0)));
                                        // System.exit(0);
                                    }
                                    partitionMSE = Double.POSITIVE_INFINITY;
                                } else {
                                    partitionMSE = obj.getMSE();
                                }

                                if (debugOptimization) {
                                    //  echoLn("\t x_" + indexSplitVariable + " Partition: " + i + " " + partitions.get(i) + " mse: " + partitionMSE);
                                }

                                //For every possible partition, we check which has the lowest MSE
                                if (partitionMSE < bestPartitionMSE || i == 0) {
                                    bestPartitionMSE = partitionMSE;
                                    optimalPartitionIndex = i;
                                    optimalZ_MSE_Right_Partition = obj.getRightMSE();
                                    optimalZ_MSE_Left_Partition = obj.getLeftMSE();
                                    numObsLeft_Partition = obj.getEffectiveNumObsLeft();
                                    numObsRight_Partition = obj.getEffectiveNumObsRight();
                                    if (debugOptimization) {
                                        echoLn("\tPartition: " + i + " " + partitions.get(i) + " mse: " + partitionMSE + " set as within-variable current best.");
                                    }
                                }
                            }

                            if (bestPartitionMSE < optimalZ_MSE || first) {
                                optimalZ = optimalPartitionIndex; // in this case, this will be the index of the best partitioning
                                // But how do we know which partition the index corresponds to?
                                optimalZ_MSE = bestPartitionMSE;
                                optimalSplitVariableIndex = indexSplitVariable;
                                optimalDiscreteCollectionIndex = collectionIndex;
                                optimalZ_MSE_Left = optimalZ_MSE_Left_Partition;
                                optimalZ_MSE_Right = optimalZ_MSE_Right_Partition;
                                numObsLeft = numObsLeft_Partition;
                                numObsRight = numObsRight_Partition;
                                first = false;
                                if (debugOptimization) {
                                    echoLn("Variable x_" + optimalSplitVariableIndex + " with partition " + partitions.get((int) optimalZ) + " giving mse of " + optimalZ_MSE + " set as overall current best.");
                                }
                            }
                        }
                    }
                }
            }

            // System.out.println("Baseline SSE is computed as: " + baseline);
            if (verbose) {
                echoLn("Number of observations in node: " + lensGrowingTree.getNumObs());
                echoLn("Improvement from " + baseline + " to " + optimalZ_MSE + " (left: " + optimalZ_MSE_Left + " [" + numObsLeft + "] right: " + optimalZ_MSE_Right + " [" + numObsRight + "])");
                echoLn("Improvement percentage (to compare against threshold): " + ((baseline - optimalZ_MSE) / baseline));
            }

            // if (baseline - optimalX_MSE < improvementThreshold) {
            if ((baseline - optimalZ_MSE) / baseline < improvementThreshold || first || baseline == 0) {
                setTerminal(true);
                if (verbose) {
                    echoLn(depth + ". Terminating due to lack of improvement in MSE; rules: " + getParentRuleDescriptive(null) + " beta: " + pmUtility.stringPrettyPrintVector(betaEstimateNode));
                }
            } else {
                setTerminal(false);
                if (discreteVector[optimalSplitVariableIndex]) {
                    ArrayList<Integer> discreteList = discreteCollection.get(optimalDiscreteCollectionIndex);
                    ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList); //Recompute all the disjoint sets of the optimal splitting variable, which is discrete

                    MomentPartitionObj obj = momentSpec.getMomentPartitionObj(lensGrowingTree, optimalSplitVariableIndex, partitions.get((int) optimalZ)); //This is where we use optimalX
                    if (verbose) {
                        echoLn(depth + ". Calculated optimal split along discrete variable, partitioning x_" + optimalSplitVariableIndex + " -> " + partitions.get((int) optimalZ) + ", generating MSE of " + obj.getMSE());
                    }
                    setRule(new SplitRule(true, optimalSplitVariableIndex, optimalZ, partitions.get((int) optimalZ)));
                    childLeft = new TreeMoment(this, momentSpec, obj.getDataSplit().getLeft(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth, null);
                    childRight = new TreeMoment(this, momentSpec, obj.getDataSplit().getRight(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth, null);
                } else {
                    MomentContinuousSplitObj obj = momentSpec.getFminObjective(lensGrowingTree, optimalSplitVariableIndex, minProportionEachPartition, minCountEachPartition);
                    if (verbose) {
                        echoLn(depth + ". Calculated optimal split along "+momentSpec.getVariableName(optimalSplitVariableIndex)+ " at " + optimalZ + ", generating MSE of " + obj.f_to_minimize(optimalZ));
                    }
                    setRule(new SplitRule(false, optimalSplitVariableIndex, optimalZ, null));

                    DataLens left = SplitContainer.getContinuousDataSplit(lensGrowingTree, optimalZ, optimalSplitVariableIndex).getLeft();
                    DataLens right = SplitContainer.getContinuousDataSplit(lensGrowingTree, optimalZ, optimalSplitVariableIndex).getRight();
                    // System.out.println("Max left: "+pmUtility.max(childLeftX, 1));
                    // System.out.println("Min right: "+pmUtility.min(childRightX, 1));
                    childLeft = new TreeMoment(this, momentSpec, left, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth, null);
                    childRight = new TreeMoment(this, momentSpec, right, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth, null);
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
         * This is the place to put some rules on what are acceptable numbers of
         * observations in order to replace the first stage estimates.
         *
         * This is also the place to put the functionality to prune the children
         * nodes if their estimates are null. That preserves that every X will
         * get an estimate.
         *
         * TODO: implement pruning if children nodes have null estimates
         */

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
}
