/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.logitRC;

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.MomentSpecification;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LogitRCDataGenerator {

    private Jama.Matrix X;
    private Jama.Matrix Y;
    private Jama.Matrix Z;

    public LogitRCDataGenerator(int numObs, MomentSpecification mySpecification, long randSeed) {
        X = new Jama.Matrix(numObs, 2);
        Z = new Jama.Matrix(numObs, mySpecification.getDiscreteVector().length);
        Y = new Jama.Matrix(numObs, 1);
        Random rng = new Random(randSeed);
        NormalDistribution normal = new NormalDistribution();
        for (int i = 0; i < numObs; i++) {
            X.set(i, 0, normal.inverse(rng.nextDouble()));
            // X.set(i, 0, 1.0);
            X.set(i, 1, Math.pow(normal.inverse(rng.nextDouble()), 2));
            // X.set(i, 1, -2.0 + 4.0*rng.nextDouble());

            // Z.set(i, 0, normal.inverse(rng.nextDouble()));
            // Z.set(i, 0, -6.0 + 10.0 * rng.nextDouble()); // uniform [-6,4]
            // Z.set(i, 1, -1.0 + 2.0 * rng.nextDouble());
//            double draw = rng.nextDouble();
//            if (draw < 0.3) {
//                Z.set(i, 2, 1);
//            } else if (draw > 0.7) {
//                Z.set(i, 2, 2);
//            }
            // Z.set(i, 0, X.get(i, 0));
            // Z.set(i, 1, X.get(i, 1));
//                pmUtility.prettyPrintVector(beta);
            Jama.Matrix subX = X.getMatrix(i, i, 0, 1); // only first two columns of X matter in producing Y
            // pmUtility.prettyPrint(subX);

            // for RC logit, want to report the aggregate share (At least to start with)
            // we will start with a dumb way; the betaTruth gets a draw from beta; we will do a stupid Monte Carlo approach to computing shares (should come back to doing quadrature if we can slot that in the framework somehow)
            boolean simpleDiscreteRC = false;
            if (simpleDiscreteRC) {
                double[][] betaList = {{-1, -1},
                {-1, 0},
                {-1, 1},
                {-1, 2},
                {-1, 3}};
                double[] weights = {0.15, 0.15, 0.3, 0.2, 0.2};
                // double[] weights = {0, 0, 1, 0, 0};
                for (int k = 0; k < betaList.length; k++) {
                    Jama.Matrix beta = new Jama.Matrix(2, 1);
                    beta.set(0, 0, betaList[k][0]);
                    beta.set(1, 0, betaList[k][1]);
                    Y.set(i, 0, Y.get(i, 0) + weights[k] * getLogitShare(subX, beta));
                }
            } else {
                double numDraws = 10000;

                // need to split on X's now
                for (int k = 0; k < numDraws; k++) {
                    Jama.Matrix beta = mySpecification.getBetaTruth(Z.getMatrix(i, i, 0, Z.getColumnDimension() - 1), rng); // Z1 and Z2 to compute beta
                    // if I don't move the beta around, I have verified that I get +Exactly+ the mean parameters back using GMM, as desired (there is no error anywhere in the DGP)
                    // adding this next line in should result in a nice RC model that should split on the variable which has a random coefficient (I hope)
                    // beta.set(0, 0, beta.get(0, 0) + normal.inverse(rng.nextDouble()));
                    beta.set(1, 0, beta.get(1, 0) + normal.inverse(rng.nextDouble()));
                    Y.set(i, 0, Y.get(i, 0) + getLogitShare(subX, beta));
                    // System.out.println(Y.get(i,0)/(k+1));
                }
                Y.set(i, 0, Y.get(i, 0) / numDraws);
            }
        }
        // pmUtility.prettyPrintVector(Y);
        // System.exit(0);
    }

    static public double getLogitShare(Jama.Matrix subX, Jama.Matrix beta) {
        double u = (subX.times(beta)).get(0, 0);
        // pmUtility.prettyPrint(subX);
        // pmUtility.prettyPrintVector(beta);
        double share = Math.exp(u) / (1.0 + Math.exp(u));
        // System.out.println("Share: "+share+" utility: "+u);
        return share;
    }

    public Matrix getX() {
        return X;
    }

    public Matrix getY() {
        return Y;
    }

    public Matrix getZ() {
        return Z;
    }

}
