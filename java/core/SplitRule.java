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

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SplitRule {

    private final boolean splitOnDiscreteVariable;
    private final int optimalSplitVariableIndex;
    private final double optimalX;
    private final IntegerPartition partition;

    /**
     *
     * This logic will encapsulate how the tree splits the sample at each node
     *
     * @param splitOnDiscreteVariable Is the split on a discrete variable?
     * @param optimalSplitVariableIndex Index of the splitting covariate
     * @param optimalX Point of X to split if continuous
     * @param partition Partition of X to split if discrete
     */
    SplitRule(boolean splitOnDiscreteVariable, int optimalSplitVariableIndex, double optimalX, IntegerPartition partition) {
        this.splitOnDiscreteVariable = splitOnDiscreteVariable;
        this.optimalSplitVariableIndex = optimalSplitVariableIndex;
        this.optimalX = optimalX;
        this.partition = partition;
    }

    @Override
    public String toString() {
        if (splitOnDiscreteVariable) {
            return partition.toString();
        }
        return "x" + optimalSplitVariableIndex + " < " + optimalX;
    }

    public int getOptimalSplitVariableIndex() {
        return optimalSplitVariableIndex;
    }

    public double getOptimalX() {
        return optimalX;
    }

    public IntegerPartition getPartition() {
        return partition;
    }

    public boolean isSplitOnDiscreteVariable() {
        return splitOnDiscreteVariable;
    }

    public boolean isLeft(Matrix xi) {
        if (splitOnDiscreteVariable) {
            Integer xc = new Integer((int) xi.get(0, optimalSplitVariableIndex));
            return partition.getLeft().contains(xc);
        } else {
            return xi.get(0, optimalSplitVariableIndex) < optimalX;
        }
    }

    String getLeftSplit() {
        if (splitOnDiscreteVariable) {
            String s = "";
            s = s.concat("x" + optimalSplitVariableIndex + " in { ");
            for (int i = 0; i < partition.getLeft().size(); i++) {
                s = s.concat(partition.getLeft().get(i) + " ");
            }
            s = s.concat("}");
            return s;
        } else {
            return "x" + optimalSplitVariableIndex + " < " + optimalX;
        }
    }

    String getLeftSplitDescriptive(MomentSpecification spec) {
        if (splitOnDiscreteVariable) {
            String s = "";
            s = s.concat(spec.getVariableName(optimalSplitVariableIndex) + " in { ");
            for (int i = 0; i < partition.getLeft().size(); i++) {
                // s = s.concat(partition.getLeft().get(i) + " ");
                s = s.concat(spec.getFixedEffectName(optimalSplitVariableIndex, partition.getLeft().get(i)) + " ");
            }
            s = s.concat("}");

            if (partition.getLeft().size() == 1) {
                s = "x"+optimalSplitVariableIndex+" = "+spec.getFixedEffectName(optimalSplitVariableIndex, partition.getLeft().get(0));
            }

            return s;
        } else {
            return spec.getVariableName(optimalSplitVariableIndex) + " < " + optimalX;
        }
    }

    SplitRuleContainer getLeftSplitContainer() {
        if (splitOnDiscreteVariable) {
            return new SplitRuleContainer(partition.getLeft(), optimalSplitVariableIndex);
        } else {
            return new SplitRuleContainer(optimalSplitVariableIndex, optimalX, true);
        }
    }

    String getRightSplit() {
        if (splitOnDiscreteVariable) {
            String s = "";
            s = s.concat("x" + optimalSplitVariableIndex + " in { ");
            for (int i = 0; i < partition.getRight().size(); i++) {
                s = s.concat(partition.getRight().get(i) + " ");
            }
            s = s.concat("}");
            return s;
        } else {
            return "x" + optimalSplitVariableIndex + " > " + optimalX;
        }
    }

    String getRightSplitDescriptive(MomentSpecification spec) {
        if (splitOnDiscreteVariable) {
            String s = "";
            s = s.concat(spec.getVariableName(optimalSplitVariableIndex) + " in { ");
            for (int i = 0; i < partition.getRight().size(); i++) {
                // s = s.concat(partition.getRight().get(i) + " ");
                s = s.concat(spec.getFixedEffectName(optimalSplitVariableIndex, partition.getRight().get(i)) + " ");
            }
            s = s.concat("}");
            
            if (partition.getRight().size() == 1) {
                s = "x"+optimalSplitVariableIndex+" = "+spec.getFixedEffectName(optimalSplitVariableIndex, partition.getRight().get(0));
            }
            
            return s;
        } else {
            return spec.getVariableName(optimalSplitVariableIndex) + " > " + optimalX;
        }
    }

    SplitRuleContainer getRightSplitContainer() {
        if (splitOnDiscreteVariable) {
            return new SplitRuleContainer(partition.getRight(), optimalSplitVariableIndex);
        } else {
            return new SplitRuleContainer(optimalSplitVariableIndex, optimalX, false);
        }
    }

}
