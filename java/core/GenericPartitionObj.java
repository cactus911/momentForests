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

/**
 * Generic partition split objective that works with any MomentSpecification.
 * Eliminates the need for model-specific partition objective subclasses by using
 * the computeOptimalBeta() factory method on MomentSpecification.
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class GenericPartitionObj extends MomentPartitionObj {

    private final SplitContainer container;
    private final MomentSpecification spec;
    private final int minCount;

    public GenericPartitionObj(IntegerPartition partition, int indexSplitVariable, DataLens lens, MomentSpecification spec, int minCount) {
        this.partition = partition;
        this.indexSplitVariable = indexSplitVariable;
        this.lens = lens;
        this.spec = spec;
        this.minCount = minCount;

        container = getDataSplit();
        numObsLeft = container.getLeft().getNumObs();
        numObsRight = container.getRight().getNumObs();
    }

    @Override
    public int getEffectiveNumObsLeft() {
        return numObsLeft;
    }

    @Override
    public int getEffectiveNumObsRight() {
        return numObsRight;
    }

    @Override
    public double getSSE() {
        leftMSE = 0;
        rightMSE = 0;

        ContainerMoment leftContainer = spec.computeOptimalBeta(container.getLeft(), false);
        ContainerMoment rightContainer = spec.computeOptimalBeta(container.getRight(), false);

        leftMSE = leftContainer.getGoodnessOfFit();
        rightMSE = rightContainer.getGoodnessOfFit();

        if (getEffectiveNumObsLeft() < minCount) {
            leftMSE = Double.POSITIVE_INFINITY;
        }
        if (getEffectiveNumObsRight() < minCount) {
            rightMSE = Double.POSITIVE_INFINITY;
        }

        return (leftMSE + rightMSE);
    }
}
