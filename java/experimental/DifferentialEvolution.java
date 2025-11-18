/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package experimental;

import JSci.maths.statistics.NormalDistribution;
import java.util.Random;
import javax.swing.ProgressMonitor;
import optimization.Uncmin_methods;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class DifferentialEvolution {

    Uncmin_methods obj;
    double optimizationTarget;

    public DifferentialEvolution(Uncmin_methods obj, double optimizationTarget) {
        this.obj = obj;
        this.optimizationTarget = optimizationTarget;
    }

    public double[] optimize(double[] guess, int numIterations, boolean reportUpdates, boolean useAdaptiveLimit) {
        int NP = 10;
        double CR = 0.9;
        double F = 0.8;

        NormalDistribution normal = new NormalDistribution();
        Random rng = new Random();

        /**
         * Initialize a cloud of starting values
         */
        double[][] currentX = new double[NP][guess.length];
        for (int i = 0; i < NP; i++) {
            for (int k = 1; k < guess.length; k++) {
                double stddev = 0.01;
                if (i == 0) {
                    stddev = 0.0;
                }
                currentX[i][k] = guess[k] + stddev * normal.inverse(rng.nextDouble());
            }
        }

        int indexBestAgent = 0;
        double fBestAgent = obj.f_to_minimize(currentX[0]);
        if (reportUpdates) {
            System.out.print(fBestAgent + " ");
            pmUtility.prettyPrint(new Jama.Matrix(currentX[0], 1));
        }
        // System.exit(0);
        boolean first = false;
        int totalIterations = 0;

        ProgressMonitor mon = new ProgressMonitor(null, "Differential Evolution", "Best: " + String.format("%.3f", fBestAgent), 1, numIterations);

        for (int t = 0; t < numIterations; t++) {
            totalIterations++;
            // loop through each agent
            for (int agent = 0; agent < NP; agent++) {
                // pick three agents at random
                int agentOne = (int) Math.floor(NP * rng.nextDouble());
                while (agent == agentOne) {
                    agentOne = (int) Math.floor(NP * rng.nextDouble());
                }
                int agentTwo = (int) Math.floor(NP * rng.nextDouble());
                while (agentTwo == agentOne || agentTwo == agent) {
                    agentTwo = (int) Math.floor(NP * rng.nextDouble());
                }
                int agentThree = (int) Math.floor(NP * rng.nextDouble());
                while (agentThree == agentOne || agentThree == agentTwo || agentThree == agent) {
                    agentThree = (int) Math.floor(NP * rng.nextDouble());
                }

                // pick a random index of the optimization problem
                int randomIndex = (int) Math.floor((guess.length - 1) * rng.nextDouble());
                // for each dimension of the optimization problem

                double[] y = new double[guess.length];

                for (int k = 1; k < guess.length; k++) {
                    // for (int k = 3; k <= 3; k++) {
                    double ri = rng.nextDouble();
                    if (ri < CR || k == randomIndex) {
                        y[k] = currentX[agentOne][k] + F * (currentX[agentTwo][k] - currentX[agentThree][k]); // evolve as a combination of other points
                    } else {
                        y[k] = currentX[agent][k]; // keep it the same
                    }
                }

                double fy = obj.f_to_minimize(y);
                double fx = obj.f_to_minimize(currentX[agent]);
                // if improvement or the same, replace current agent's guess
                if (fy <= fx) {
                    for (int k = 1; k < guess.length; k++) {
                        currentX[agent][k] = y[k];
                    }
                }

                if (fy <= fBestAgent || first) {
                    double pctImprovement = (fBestAgent - fy) / fBestAgent;
                    
                    fBestAgent = fy;
                    indexBestAgent = agent;
                    first = false;

                    if (reportUpdates) {
                        System.out.format("[%d, %.2f%%] %f ", totalIterations, pctImprovement, fy);
                        pmUtility.prettyPrint(new Jama.Matrix(y, 1));
                    }

                    // don't let this stop if it is still making progress
                    if (useAdaptiveLimit && pctImprovement>1E-2) {
                        // t = Math.min(t, numIterations - 50);
                        t = Math.min(t, numIterations/2);
                    }
                }

                if (fBestAgent <= optimizationTarget) {
                    agent = NP;
                    t = numIterations;
                }
            }
            mon.setProgress(t + 1);
            // mon.setNote((t + 1) + " / " + numIterations);
            mon.setNote("Best: " + String.format("%.3f", fBestAgent));
        }

        return currentX[indexBestAgent];
    }

}
