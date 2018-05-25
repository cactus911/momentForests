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
import optimization.Fmin;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class TreeMoment {

    TreeMoment parent;
    Jama.Matrix nodeX;
    Jama.Matrix nodeY;
    TreeMoment childLeft;
    TreeMoment childRight;
    Boolean isLeftNode;
    Boolean terminal = true;
    Boolean[] discreteVector;

    SplitRule rule;
    private int depth;

    MomentSpecification momentSpec;
    Jama.Matrix betaEstimateNode;
    Jama.Matrix varianceMatrix;

    private final ArrayList<Jama.Matrix> honestXList = new ArrayList<>();
    private final ArrayList<Jama.Matrix> honestYList = new ArrayList<>();
    private Jama.Matrix honestX = null;
    private Jama.Matrix honestY = null;
    int numHonestXObservations;

    // private final static double improvementThreshold = 0.01;
    int maxDepth = 100;
    double minProportionEachPartition;
    int minCountEachPartition;
    double improvementThreshold;

    boolean verbose;

    boolean debugOptimization = false;

    public TreeMoment(TreeMoment parent, MomentSpecification spec, Matrix subsampleX, Matrix subsampleY,
            Boolean[] discreteVector, boolean verbose, double minProportionEachPartition,
            int minCountEachPartition, double improvementThreshold, boolean isLeft, int maxDepth) {
        this.momentSpec = spec;
        this.parent = parent;
        this.nodeX = subsampleX;
        this.nodeY = subsampleY;
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
     * Return tree's prediction of treatment effect on the basis of x_i.
     *
     * @param xi
     * @return
     */
    public Double getPredictedY(Jama.Matrix xi) {
        if (terminal) {
            // return the treatment effect for this x_i, not the predicted y_i
            return momentSpec.getPredictedY(xi, betaEstimateNode);
        } else if (rule.isLeft(xi)) {
            return childLeft.getPredictedY(xi);
        } else {
            return childRight.getPredictedY(xi);
        }
    }

    public void setDepth(int depth) {
        this.depth = depth;
    }

    public int getDepth() {
        return depth;
    }

    public Matrix getSubsampleX() {
        return nodeX;
    }

    public void setSubsampleY(Matrix subsampleY) {
        this.nodeY = subsampleY;
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

    public void determineSplit() {
        if (verbose) {
            System.out.println("----------- computing optimal split --------------");
            System.out.println(depth + " " + getParentRule(null));
        }

//        System.out.println("in this sub-node:");
//        if (nodeY.getRowDimension() > 10) {
//            pmUtility.prettyPrint(pmUtility.concatMatrix(nodeY, nodeX).getMatrix(0, 10, 0, nodeX.getColumnDimension()));
//            System.out.println("TreeMoment.java:152 -> min: "+pmUtility.min(nodeX, 1)+" max: "+pmUtility.max(nodeX, 1));
//        }
        /**
         * This is the place to set a priori conditions on growing the tree (max
         * depth, etc.)
         */
        if (depth < maxDepth) {
            double optimalX = 0;
            double optimalX_MSE = 0;
            double optimalX_MSE_Left = 0;
            double optimalX_MSE_Right = 0;
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
            for (int k : momentSpec.getVariableIndicesToSearchOver()) {
                if (discreteVector[k] == true) {
                    TreeSet<Integer> discreteTreeSet = new TreeSet<>();
                    for (int i = 0; i < nodeX.getRowDimension(); i++) {
                        discreteTreeSet.add((int) nodeX.get(i, k));
                    }
                    ArrayList<Integer> discreteList = new ArrayList<>(discreteTreeSet);
                    discreteCollection.add(discreteList);
                    // System.out.println("Added list for k = "+k+" of size "+discreteList.size());
                    discreteCollectionIndex.add(k);
                }
            }

            /**
             * We are going to store p integers here that we will search over.
             * This is going to be random at each split, we don't have to keep
             * track across branches.
             */
            TreeSet<Integer> randomForestIndex = new TreeSet<>();
            boolean useRandomForest = false;
            if (useRandomForest) {
                int P = Math.min(100, momentSpec.getVariableIndicesToSearchOver().length);
                for (int i = 0; i < P; i++) {
                    int index = (int) Math.floor(Math.random() * momentSpec.getVariableIndicesToSearchOver().length);
                    while (randomForestIndex.contains(index)) {
                        index = (int) Math.floor(Math.random() * momentSpec.getVariableIndicesToSearchOver().length);
                    }
                    randomForestIndex.add(index);
                }
            } else {
                for (int index: momentSpec.getVariableIndicesToSearchOver()) {
                    randomForestIndex.add(index);
                }
            }

            boolean first = true;

            for (int indexSplitVariable : momentSpec.getVariableIndicesToSearchOver()) {
                if(debugOptimization) {
                    System.out.println("indexSplitVariable: "+indexSplitVariable+" isDiscrete: "+discreteVector[indexSplitVariable]+" "+randomForestIndex.contains(indexSplitVariable));
                }
                if (randomForestIndex.contains(indexSplitVariable)) {
                    if (discreteVector[indexSplitVariable] == false) {
                        MomentContinuousSplitObj obj = momentSpec.getFminObjective(nodeY, nodeX, indexSplitVariable, minProportionEachPartition, minCountEachPartition);
                        double minX = pmUtility.min(nodeX, indexSplitVariable);
                        double maxX = pmUtility.max(nodeX, indexSplitVariable);
                        // System.out.println("TreeMoment.java:224 -> min x_1: "+minX+" max x_1: "+maxX);
                        double optimalX_k = Fmin.fmin(minX, maxX, obj, 1E-16);
                        double optimalX_MSE_k = obj.f_to_minimize(optimalX_k);

                        if (debugOptimization) {
                            System.out.println("\tFmin search on x_" + indexSplitVariable + " found x = " + optimalX_k + " mse: " + optimalX_MSE_k);
                        }

                        boolean testGridSearch = true;
                        double h = 1E-5;
                        if (testGridSearch) {
                            for (double x = minX; x <= maxX; x += h + (maxX - minX) / 100.0) {
                                double f = obj.f_to_minimize(x);
                                if (debugOptimization) {
                                    System.out.println("\tGrid search x_" + indexSplitVariable + " from " + optimalX_MSE_k + " to " + f + " by moving from " + optimalX_k + " to " + x);
                                }
                                if (f < optimalX_MSE_k) {
                                    optimalX_k = x;
                                    optimalX_MSE_k = obj.f_to_minimize(x);
                                }
                            }
                        }

                        if (optimalX_MSE_k < optimalX_MSE || first) {
                            optimalX = optimalX_k;
                            optimalX_MSE = optimalX_MSE_k;
                            optimalSplitVariableIndex = indexSplitVariable;
                            optimalX_MSE_Left = obj.getLeftMSE();
                            optimalX_MSE_Right = obj.getRightMSE();
                            numObsLeft = obj.getNumObsLeft();
                            numObsRight = obj.getNumObsRight();
                            first = false;
                        }
                    } else {
                        int collectionIndex = discreteCollectionIndex.indexOf(indexSplitVariable);
                        ArrayList<Integer> discreteList = discreteCollection.get(collectionIndex);
                        ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList);
                        // System.out.println("Partition size: "+partitions.size()+" discreteList.size(): "+discreteList.size());
                        /**
                         * Need to put in a check here that the
                         * discreteList.size is greater than one element
                         */
                        if (discreteList.size() > 1) {
                            int optimalPartitionIndex = 0;
                            double bestPartitionMSE = 0;
                            double optimalX_MSE_Right_Partition = 0;
                            double optimalX_MSE_Left_Partition = 0;
                            int numObsLeft_Partition = 0;
                            int numObsRight_Partition = 0;
                            for (int i = 0; i < partitions.size(); i++) {
                                MomentPartitionObj obj = momentSpec.getMomentPartitionObj(nodeX, nodeY, indexSplitVariable, partitions.get(i));

                                double partitionMSE = 0;
                                if (obj.getNumObsLeft() < minCountEachPartition || obj.getNumObsRight() < minCountEachPartition) {
                                    if (debugOptimization) {
                                        System.out.println("\t\tMin K violated: rejecting partition for left obs: " + obj.getNumObsLeft() + " right obs: " + obj.getNumObsRight());
                                    }
                                    partitionMSE = Double.POSITIVE_INFINITY;
                                } else if (((obj.getNumObsLeft() + 0.0) / (nodeX.getRowDimension() + 0.0)) < minProportionEachPartition || ((obj.getNumObsRight() + 0.0) / (nodeX.getRowDimension() + 0.0)) < minProportionEachPartition) {
                                    if (debugOptimization) {
                                        System.out.println("\t\tRejecting partition for proportion; left: " + ((obj.getNumObsLeft() + 0.0) / (nodeX.getRowDimension() + 0.0)) + " right: " + ((obj.getNumObsRight() + 0.0) / (nodeX.getRowDimension() + 0.0)));
                                        // System.exit(0);
                                    }
                                    partitionMSE = Double.POSITIVE_INFINITY;
                                } else {
                                    partitionMSE = obj.getMSE();
                                }

                                if (debugOptimization) {
                                    System.out.println("\t x_" + indexSplitVariable + " Partition: " + i + " " + partitions.get(i) + " mse: " + partitionMSE);
                                }

                                if (partitionMSE < bestPartitionMSE || i == 0) {
                                    bestPartitionMSE = partitionMSE;
                                    optimalPartitionIndex = i;
                                    optimalX_MSE_Right_Partition = obj.getRightMSE();
                                    optimalX_MSE_Left_Partition = obj.getLeftMSE();
                                    numObsLeft_Partition = obj.getNumObsLeft();
                                    numObsRight_Partition = obj.getNumObsRight();
                                    if (debugOptimization) {
                                        System.out.println("\tPartition: " + i + " " + partitions.get(i) + " mse: " + partitionMSE + " set as within-variable current best.");
                                    }
                                }
                            }
                            if (bestPartitionMSE < optimalX_MSE || first) {
                                optimalX = optimalPartitionIndex; // in this case, this will be the index of the best partitioning
                                optimalX_MSE = bestPartitionMSE;
                                optimalSplitVariableIndex = indexSplitVariable;
                                optimalDiscreteCollectionIndex = collectionIndex;
                                optimalX_MSE_Left = optimalX_MSE_Left_Partition;
                                optimalX_MSE_Right = optimalX_MSE_Right_Partition;
                                numObsLeft = numObsLeft_Partition;
                                numObsRight = numObsRight_Partition;
                                first = false;
                                if (debugOptimization) {
                                    System.out.println("Variable x_" + optimalSplitVariableIndex + " with partition " + partitions.get((int) optimalX) + " giving mse of " + optimalX_MSE + " set as overall current best.");
                                }
                            }
                        }
                    }
                }
            }

            /**
             * Now compute a baseline MSE for comparing the improvement in fit
             */
            // System.out.println("Computing baseline SSE");
            ContainerMoment currentNodeMoment = momentSpec.computeOptimalBeta(nodeY, nodeX);
            setNodeEstimatedBeta(currentNodeMoment.getBeta());
            setNodeEstimatedVariance(currentNodeMoment.getVariance());

            double baseline = currentNodeMoment.getMSE();
            // System.out.println("Baseline SSE is computed as: " + baseline);
            if (verbose) {
                System.out.println("Improvement from " + baseline + " to " + optimalX_MSE + " (left: " + optimalX_MSE_Left + " [" + numObsLeft + "] right: " + optimalX_MSE_Right + " [" + numObsRight + "])");
            }

            /**
             * There may be a problem here in that first is never hit when the
             * first index of X is not evaluated. The basic issue is that the
             * index of X=0 is set as the default best. That leads to insanity
             * when you are doing discrete cuts on X1 and ignoring X0. Need to
             * add in check here for whether first is still true (which means
             * that there were no valid splits, and we should terminate)
             */
            // if (baseline - optimalX_MSE < improvementThreshold) {
            if ((baseline - optimalX_MSE) / baseline < improvementThreshold || first || baseline == 0) {
                setTerminal(true);
                if (verbose) {
                    System.out.println(depth + ". Terminating due to lack of improvement in MSE; rules: " + getParentRule(null) + " beta: " + pmUtility.stringPrettyPrintVector(betaEstimateNode));
                }
            } else {
                setTerminal(false);
                if (discreteVector[optimalSplitVariableIndex]) {
                    ArrayList<Integer> discreteList = discreteCollection.get(optimalDiscreteCollectionIndex);
                    ArrayList<IntegerPartition> partitions = DisjointSet.computeAllDisjointSets(discreteList);

                    MomentPartitionObj obj = momentSpec.getMomentPartitionObj(nodeX, nodeY, optimalSplitVariableIndex, partitions.get((int) optimalX));
                    if (verbose) {
                        System.out.println(depth + ". Calculated optimal split along discrete variable, partitioning x_" + optimalSplitVariableIndex + " -> " + partitions.get((int) optimalX) + ", generating MSE of " + obj.getMSE());
                    }
                    rule = new SplitRule(true, optimalSplitVariableIndex, optimalX, partitions.get((int) optimalX));
                    childLeft = new TreeMoment(this, momentSpec, obj.getDataSplit().getxLeft(), obj.getDataSplit().getyLeft(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth);
                    childRight = new TreeMoment(this, momentSpec, obj.getDataSplit().getxRight(), obj.getDataSplit().getyRight(), discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth);
                } else {
                    MomentContinuousSplitObj obj = momentSpec.getFminObjective(nodeY, nodeX, optimalSplitVariableIndex, minProportionEachPartition, minCountEachPartition);
                    if (verbose) {
                        System.out.println(depth + ". Calculated optimal split along variable " + optimalSplitVariableIndex + " at " + optimalX + ", generating MSE of " + obj.f_to_minimize(optimalX));
                    }
                    rule = new SplitRule(false, optimalSplitVariableIndex, optimalX, null);

                    Jama.Matrix childLeftX = SplitContainer.getContinuousDataSplit(nodeY, nodeX, optimalX, optimalSplitVariableIndex).getxLeft();
                    Jama.Matrix childRightX = SplitContainer.getContinuousDataSplit(nodeY, nodeX, optimalX, optimalSplitVariableIndex).getxRight();
                    Jama.Matrix childLeftY = SplitContainer.getContinuousDataSplit(nodeY, nodeX, optimalX, optimalSplitVariableIndex).getyLeft();
                    Jama.Matrix childRightY = SplitContainer.getContinuousDataSplit(nodeY, nodeX, optimalX, optimalSplitVariableIndex).getyRight();
                    // System.out.println("Max left: "+pmUtility.max(childLeftX, 1));
                    // System.out.println("Min right: "+pmUtility.min(childRightX, 1));
                    childLeft = new TreeMoment(this, momentSpec, childLeftX, childLeftY, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            true, maxDepth);
                    childRight = new TreeMoment(this, momentSpec, childRightX, childRightY, discreteVector, verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold,
                            false, maxDepth);
                }
                childLeft.determineSplit();
                childRight.determineSplit();
            }
        } else {
            setTerminal(true);
            // System.out.print("before call: ");
            // pmUtility.prettyPrintVector(betaEstimateNode);
            ContainerMoment currentNodeMoment = momentSpec.computeOptimalBeta(nodeY, nodeX);
            setNodeEstimatedBeta(currentNodeMoment.getBeta());
            // System.out.print("after call: ");
            // pmUtility.prettyPrintVector(betaEstimateNode);
            setNodeEstimatedVariance(currentNodeMoment.getVariance());

            if (verbose) {
                System.out.println(depth + ". Terminal beta: " + pmUtility.stringPrettyPrintVector(betaEstimateNode));
                // System.out.println(depth + ". Terminal RDD value: " + getRDDEstimate());
            }
        }
    }

    public String getParentRule(TreeSet<Integer> indexPreviousSplits) {
        if (terminal && parent == null) {
            return "{ Stump }";
        }
        if (parent == null) {
            return " }";
        }

        SplitRule r = parent.rule;
        String s = "";
        if (terminal) {
            s = "{ ";
        }
        if (r.isSplitOnDiscreteVariable()) {
            if (isLeftNode) {
                s = s + r.getLeftSplit();
            } else {
                s = s + r.getRightSplit();
            }
        } else if (isLeftNode) {
            s = s + "x_" + r.getOptimalSplitVariableIndex() + " < " + r.getOptimalX();
        } else {
            s = s + "x_" + r.getOptimalSplitVariableIndex() + " > " + r.getOptimalX();
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
        return s + parent.getParentRule(indexPreviousSplits);
    }

    public String getParentRuleDescriptive(TreeSet<Integer> indexPreviousSplits) {
        if (terminal && parent == null) {
            return "{ Stump }";
        }
        if (parent == null) {
            return " }";
        }

        SplitRule r = parent.rule;
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
            s = s + momentSpec.getVariableName(r.getOptimalSplitVariableIndex()) + " < " + r.getOptimalX();
        } else {
            s = s + momentSpec.getVariableName(r.getOptimalSplitVariableIndex()) + " > " + r.getOptimalX();
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

    public void setNodeEstimatedBeta(Matrix beta) {
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
            ruleList.add(rule);
            childLeft.addAllRulesRecursively(ruleList);
            childRight.addAllRulesRecursively(ruleList);
        }
    }

    public void sortXToCorrectLeafs(Matrix yi, Matrix xi) {
        honestXList.add(xi.copy());
        honestYList.add(yi.copy());
        if (!terminal) {
            if (rule.isLeft(xi)) {
                // System.out.println("Went left");
                childLeft.sortXToCorrectLeafs(yi, xi);
            } else {
                // System.out.println("Went right");
                childRight.sortXToCorrectLeafs(yi, xi);
            }
        }
    }

    public void consolidateHonestData() {
        if (honestXList.size() > 0) {
            honestX = new Jama.Matrix(honestXList.size(), honestXList.get(0).getColumnDimension());
            honestY = new Jama.Matrix(honestXList.size(), 1);
            for (int i = 0; i < honestXList.size(); i++) {
                for (int j = 0; j < honestXList.get(0).getColumnDimension(); j++) {
                    honestX.set(i, j, honestXList.get(i).get(0, j));
                    honestY.set(i, 0, honestYList.get(i).get(0, 0));
                }
            }
        }
        if (!terminal) {
            childLeft.consolidateHonestData();
            childRight.consolidateHonestData();
        }
    }

    public void printTree() {
        if (!terminal) {
            childLeft.printTree();
            childRight.printTree();
        } else {
            // System.out.println(getParentRuleDescriptive(null) + " " + pmUtility.stringPrettyPrintVector(betaEstimateNode) + " (" + pmUtility.stringPrettyPrintVector(varianceMatrix) + ")");
            System.out.println(getParentRuleDescriptive(null) + " ["+numHonestXObservations+"] " + momentSpec.formatTreeLeafOutput(betaEstimateNode, varianceMatrix));
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
        if (terminal) {
            if (verbose) {
                System.out.print(getParentRule(null) + " honest set to ");
            }
            if (honestY == null) {
                setNodeEstimatedBeta(null);
                setNodeEstimatedVariance(null);
                if (verbose) {
                    System.out.println("honestY null");
                }
            } else {
                ContainerMoment c = momentSpec.computeOptimalBeta(honestY, honestX);
                setNodeEstimatedBeta(c.getBeta());
                setNodeEstimatedVariance(c.getVariance());
                if (verbose) {
//                    System.out.println("beta:");
//                    pmUtility.prettyPrintVector(c.getBeta());
//                    System.out.println("var(beta):");
//                    pmUtility.prettyPrint(c.getVariance());
                    // System.exit(0);

                    if (c.getBeta() != null) {
                        System.out.println(c.getBeta().get(0, 0) + " (" + Math.sqrt(c.getVariance().get(0, 0)) + ")");
                    } else {
                        System.out.println("null");
                    }
                }
            }
        } else {
            childLeft.estimateHonestTree();
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
                            System.out.println("Detected child left is null");
                        }
                    }
                    if (rightBeta == null) {
                        if (verbose) {
                            System.out.println("Detected child right is null");
                        }
                    }
                    setTerminal(true);
                    if (honestX == null) {
                        setNodeEstimatedBeta(null);
                        setNodeEstimatedVariance(null);
                        // System.out.println("null pruned");
                    } else {
                        ContainerMoment c = momentSpec.computeOptimalBeta(honestY, honestX);
                        setNodeEstimatedBeta(c.getBeta());
                        setNodeEstimatedVariance(c.getVariance());
                        if (verbose) {
                            System.out.println("pruned children, replaced with: n = " + honestX.getRowDimension() + " " + pmUtility.stringPrettyPrintVector(c.getBeta()));
                        }
                    }
                } else {
                    // System.out.println("not pruned: n = " + honestX.getRowDimension());
                }
            }
        }
    }

    public void clearEstimationData() {
        nodeX = null;
        nodeY = null;
        if (!terminal) {
            childLeft.clearEstimationData();
            childRight.clearEstimationData();
        }
    }

    public void clearHonestyData() {
        numHonestXObservations = honestX.getRowDimension();
        honestX = null;
        honestY = null;
        if (!terminal) {
            childLeft.clearHonestyData();
            childRight.clearHonestyData();
        }
    }

    public Matrix getEstimatedBeta(Matrix xi) {
        if (terminal) {
            // return the treatment effect for this x_i, not the predicted y_i
            return betaEstimateNode;
        } else // use rule to figure out whether to return left or right node's value
        // this will keep going down the rabbit hole until it returns a terminal node's value
        // kind of cool how this works!
         if (rule.isLeft(xi)) {
                return childLeft.getEstimatedBeta(xi);
            } else {
                return childRight.getEstimatedBeta(xi);
            }
    }

    private void setNodeEstimatedVariance(Matrix variance) {
        varianceMatrix = variance;
    }

    public Matrix getVarianceMatrix(Jama.Matrix xi) {
        if (terminal) {
            // return the treatment effect for this x_i, not the predicted y_i
            return varianceMatrix;
        } else // use rule to figure out whether to return left or right node's value
        // this will keep going down the rabbit hole until it returns a terminal node's value
        // kind of cool how this works!
         if (rule.isLeft(xi)) {
                return childLeft.getVarianceMatrix(xi);
            } else {
                return childRight.getVarianceMatrix(xi);
            }
    }

}
