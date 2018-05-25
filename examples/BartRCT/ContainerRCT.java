/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.BartRCT;

import Jama.Matrix;
import core.ContainerMoment;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class ContainerRCT extends ContainerMoment {

    double mse;
    Jama.Matrix beta;
    Jama.Matrix variance;

    public ContainerRCT(Jama.Matrix X, Jama.Matrix Y) {
        computeBetaAndErrors(X, Y);
    }

    private void computeBetaAndErrors(Jama.Matrix X, Jama.Matrix Y) {
        double meanTreatment = 0;
        double meanControl = 0;
        double varianceTreatment = 0;
        double varianceControl = 0;

        int countTreatment = 0;
        int countControl = 0;
        for (int i = 0; i < Y.getRowDimension(); i++) {
            if (X.get(i, 0) == 0) {
                meanControl += Y.get(i, 0);
                countControl++;
            } else {
                meanTreatment += Y.get(i, 0);
                countTreatment++;
            }
        }

        meanControl /= countControl;
        meanTreatment /= countTreatment;

        for (int i = 0; i < Y.getRowDimension(); i++) {
            if (X.get(i, 0) == 0) {
                varianceControl += Math.pow(Y.get(i, 0) - meanControl, 2);
            } else {
                varianceTreatment += Math.pow(Y.get(i, 0) - meanTreatment, 2);
            }
        }
        varianceControl /= countControl;
        varianceTreatment /= countTreatment;
        // System.out.format("Treatment Mean: %f Variance: %f %n", meanTreatment, varianceTreatment);
        // System.out.format("  Control Mean: %f Variance: %f %n", meanControl, varianceControl);

        if (countControl == 0 || countTreatment == 0) {
            // System.out.println("Setting beta to NULL");
            // System.out.println("numObs_control: "+countControl+" numObs_treatment: "+countTreatment);
            beta = null;
            variance = null;
            mse = Double.POSITIVE_INFINITY;
        } else {
            beta = new Jama.Matrix(1, 1, meanTreatment - meanControl);
            variance = new Jama.Matrix(1, 1, (varianceControl / countControl) + (varianceTreatment / countTreatment));

            // pmUtility.prettyPrintVector(Y);
            double sse = 0;
            for (int i = 0; i < Y.getRowDimension(); i++) {
                if (X.get(i, 0) == 0) {
                    sse += Math.pow(Y.get(i, 0) - meanControl, 2);
                } else {
                    sse += Math.pow(Y.get(i, 0) - meanTreatment, 2);
                }
            }
            mse = sse; // / Y.getRowDimension();
            // System.out.println("containerRCT.java:77, mse = "+mse);

            /**
             * Try this with a simple OLS and compare just to make sure
             *
             * I verified that these give identical answers.
             */
            boolean tryOLS = false;
            if (tryOLS) {
                System.out.println(" n_control: " + countControl + " n_treatment: " + countTreatment);
                Jama.Matrix W = pmUtility.getColumn(X, 0);
                Jama.Matrix betaOLS = pmUtility.OLS(W, Y, true);
                pmUtility.prettyPrint(pmUtility.concatMatrix(Y, W));
                System.out.print("beta_hand: " + beta.get(0, 0) + " ");
                pmUtility.prettyPrintVector(betaOLS);
                double sseOLS = 0;
                for (int i = 0; i < Y.getRowDimension(); i++) {
                    double fitY = betaOLS.get(0, 0) + betaOLS.get(1, 0) * W.get(i, 0);
                    sseOLS += Math.pow(Y.get(i, 0) - fitY, 2);
                }
                System.out.println("sse: " + sse + " sseOLS: " + sseOLS);
            }
            // System.exit(0);
        }
    }

    @Override
    public Matrix getBeta() {
        return beta;
    }

    @Override
    public double getMSE() {
        return mse;
    }

    @Override
    public Matrix getVariance() {
        return variance;
    }

}
