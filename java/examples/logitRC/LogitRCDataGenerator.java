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
            // X.set(i, 0, normal.inverse(rng.nextDouble()));
            X.set(i, 0, 1.0);
            X.set(i, 1, Math.pow(normal.inverse(rng.nextDouble()), 2));
            // X.set(i, 1, -2.0 + 4.0*rng.nextDouble());

            Z.set(i, 0, normal.inverse(rng.nextDouble()));
            // Z.set(i, 0, -6.0 + 10.0 * rng.nextDouble()); // uniform [-6,4]
            Z.set(i, 1, -1.0 + 2.0 * rng.nextDouble());

            double draw = rng.nextDouble();
            if (draw < 0.3) {
                Z.set(i, 2, 1);
            } else if (draw > 0.7) {
                Z.set(i, 2, 2);
            }
            // Z.set(i, 0, X.get(i, 0));
            // Z.set(i, 1, X.get(i, 1));

//                pmUtility.prettyPrintVector(beta);
            Jama.Matrix subX = X.getMatrix(i, i, 0, 1); // only first two columns of X matter in producing Y
            // pmUtility.prettyPrint(subX);

            // for RC logit, want to report the aggregate share (At least to start with)
            // we will start with a dumb way; the betaTruth gets a draw from beta; we will do a stupid Monte Carlo approach to computing shares (should come back to doing quadrature if we can slot that in the framework somehow)
            Jama.Matrix beta = mySpecification.getBetaTruth(Z.getMatrix(i, i, 0, Z.getColumnDimension() - 1), rng); // Z1 and Z2 to compute beta
            Y.set(i, 0, getLogitShare(subX, beta, rng));
        }
    }

    static public double getLogitShare(Jama.Matrix subX, Jama.Matrix beta, Random rng) {
        double u = (subX.times(beta)).get(0, 0);
        double share = Math.exp(u) / (1.0 + Math.exp(u));
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
