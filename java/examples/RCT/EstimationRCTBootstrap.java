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
package examples.RCT;

import com.stata.sfi.*;
import JSci.maths.statistics.NormalDistribution;
import core.TreeMoment;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;
import java.util.TreeSet;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class EstimationRCTBootstrap {
	
	 public static int Momentforests(String[] args) {  
             
		 double resvv;
		 long bootsobs = Data.getObsParsedIn2();
		 resvv = Data.getNum(1, bootsobs);  
		 int resv = (int) resvv ;
        
                int rc = 0;
		EstimationRCTBootstrap go = new EstimationRCTBootstrap(resv);
		return(rc);
     }
    

         
       @SuppressWarnings("deprecation")
	public EstimationRCTBootstrap(int resv) {
   		MomentSpecificationRCTbootstrap bartbootstrap = new MomentSpecificationRCTbootstrap();
	    bartbootstrap.loadData();		 
        

        int bestK = 0;
        double bestAlpha = 0;
        double bestMSEBar = 0;
        int bestDepth = 0;
        boolean first = true;
        double bestMSPE = 0;
        double bestProportion = 0;
                
        int numObs = bartbootstrap.getX().getRowDimension();
        Jama.Matrix allX = bartbootstrap.getX().copy();
        Jama.Matrix allXoriginal = bartbootstrap.getXoriginal().copy();

        // try implementing a true OOS prediction error here?
        for (double proportion = 0.3; proportion <= 0.3; proportion += 0.05) {
            int numObsEstimateLeafValues = (int) Math.floor(numObs * proportion);
            int numPredictObs = (int) Math.floor(numObs * 0.1); // it was originally 10
            int numObsGrowTreeStructure = numObs - numObsEstimateLeafValues - numPredictObs; // halfObs / 2;
            SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsEstimateLeafValues + " predictObs: " + numPredictObs);

            Jama.Matrix treeX = new Jama.Matrix(numObsGrowTreeStructure, bartbootstrap.getX().getColumnDimension());
            Jama.Matrix treeY = new Jama.Matrix(numObsGrowTreeStructure, 1);
            Jama.Matrix honestX = new Jama.Matrix(numObsEstimateLeafValues, bartbootstrap.getX().getColumnDimension());
            Jama.Matrix honestY = new Jama.Matrix(numObsEstimateLeafValues, 1);
            Jama.Matrix predictX = new Jama.Matrix(numPredictObs, bartbootstrap.getX().getColumnDimension());
            Jama.Matrix predictY = new Jama.Matrix(numPredictObs, 1);

            TreeSet<Integer> growTreeSet = new TreeSet<>();
            int count = 0;
            while (count < numObsGrowTreeStructure) {
                if (count < numObsGrowTreeStructure) {
                    int index = (int) Math.floor(Math.random() * numObs);
                    if (!growTreeSet.contains(index)) {
                        growTreeSet.add(index);
                        count++;
                    }
                }
            }
            TreeSet<Integer> predictSet = new TreeSet<>();
            count = 0;
            while (count < numPredictObs) {
                if (count < numPredictObs) {
                    int index = (int) Math.floor(Math.random() * numObs);
                    if (!predictSet.contains(index) && !growTreeSet.contains(index)) {
                        predictSet.add(index);
                        count++;
                    }
                }
            }
            SFIToolkit.displayln("honestX row: " + honestX.getRowDimension() + " col: " + honestX.getColumnDimension());
            SFIToolkit.displayln("treeSet.size(): " + growTreeSet.size());

            int countTree = 0;
            int countPredict = 0;
            int countHonest = 0;
            for (int i = 0; i < numObs; i++) {
                if (growTreeSet.contains(i)) {
                    for (int j = 0; j < bartbootstrap.getX().getColumnDimension(); j++) {
                        treeX.set(countTree, j, bartbootstrap.getX().get(i, j));
                    }
                    treeY.set(countTree, 0, bartbootstrap.getY().get(i, 0));  
                    countTree++;
                } else if (predictSet.contains(i)) {
                    for (int j = 0; j < bartbootstrap.getX().getColumnDimension(); j++) {
                        predictX.set(countPredict, j, bartbootstrap.getX().get(i, j));
                    }
                    predictY.set(countPredict, 0, bartbootstrap.getY().get(i, 0));
                    countPredict++;
                } else {
                    for (int j = 0; j < bartbootstrap.getX().getColumnDimension(); j++) {
                        honestX.set(countHonest, j, bartbootstrap.getX().get(i, j));
                    }
                    honestY.set(countHonest, 0, bartbootstrap.getY().get(i, 0));
                    countHonest++;
                }
                // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
            }

            double minProportionEachPartition = 0.001;
            int minCountEachPartition = 10;
            double improvementThreshold = 0.01;
            int maxDepth = 100;
            boolean verbose = false;
            
            Jama.Matrix CVParameters = bartbootstrap.cvparameters();
               

            // for (maxDepth = 1; maxDepth <= 20; maxDepth++) {
            for (improvementThreshold = CVParameters.get(0, 1); improvementThreshold <= CVParameters.get(0, 1); improvementThreshold *= 100) {
                for (minCountEachPartition = (int) CVParameters.get(0, 0); minCountEachPartition <= (int) CVParameters.get(0, 0); minCountEachPartition += 50) {
                	SFIToolkit.displayln("MSE_bar: " + improvementThreshold + " k: " + minCountEachPartition);
                    TreeMoment momentTree = new TreeMoment(null, bartbootstrap, treeX, treeY, bartbootstrap.getDiscreteVector(), verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth);
                    momentTree.determineSplit();
                    momentTree.printTree();

                    double MSE = 0;
                    double MSPE = 0;
                    double counter = 0;
                    int nullCounter = 0;
                    for (int i = 0; i < treeX.getRowDimension(); i++) {
                        Double predictedY = momentTree.getPredictedY(treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1));
                        if (predictedY != null) {
                            MSE += Math.pow(treeY.get(i, 0) - predictedY, 2);
                            counter++;
                        } else {
                            pmUtility.prettyPrint(treeX.getMatrix(i, i, 0, treeX.getColumnDimension() - 1));
                            nullCounter++;
                        }
                    }
                    MSE /= counter;

                    boolean computeHonestTree = true;
                    if (computeHonestTree) {
                        boolean TIMING_DEBUG = false;

                        long t1 = System.currentTimeMillis();
                        for (int i = 0; i < honestX.getRowDimension(); i++) {
                            Jama.Matrix xi = honestX.getMatrix(i, i, 0, honestX.getColumnDimension() - 1);
                            Jama.Matrix yi = honestY.getMatrix(i, i, 0, honestY.getColumnDimension() - 1);
                            momentTree.sortXToCorrectLeafs(yi, xi);
                        }
                        momentTree.consolidateHonestData();
                        long t2 = System.currentTimeMillis();
                        if (TIMING_DEBUG) {
                        	SFIToolkit.displayln("Sorted observations into tree in " + (t2 - t1) + " ms.");
                        }

                        t1 = System.currentTimeMillis();
                        momentTree.estimateHonestTree();
                        t2 = System.currentTimeMillis();
                        if (TIMING_DEBUG) {
                        	SFIToolkit.displayln("Recomputed estimates in " + (t2 - t1) + " ms.");
                        }
                    }

                    counter = 0;
                    int nullCounterMSPE = 0;
                    for (int i = 0; i < predictX.getRowDimension(); i++) {
                        Double predictedY = momentTree.getPredictedY(predictX.getMatrix(i, i, 0, predictX.getColumnDimension() - 1));
                        if (predictedY != null) {
                            MSPE += Math.pow(predictY.get(i, 0) - predictedY, 2);
                            counter++;
                        } else {
                            pmUtility.prettyPrint(predictX.getMatrix(i, i, 0, honestX.getColumnDimension() - 1));
                            nullCounterMSPE++;
                        }
                    }
                    MSPE /= counter;

                    SFIToolkit.displayln("maxDepth: " + maxDepth + " " + " k: " + minCountEachPartition + " " + " alpha: " + minProportionEachPartition + " " + " mse_bar: " + improvementThreshold + " ");
                    SFIToolkit.displayln("MSE: " + MSE + " MSPE: " + MSPE + " nulls: " + nullCounter + " ");
                    // System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);
                    if (MSPE < bestMSPE || first) {
                    	SFIToolkit.displayln(" * ");
                        bestK = minCountEachPartition;
                        bestAlpha = minProportionEachPartition;
                        bestMSEBar = improvementThreshold;
                        bestDepth = maxDepth;
                        bestMSPE = MSPE;
                        bestProportion = proportion;
                        first = false;
                    } else {
                    	SFIToolkit.displayln(" ");
                    }
                }
            }
        }
        SFIToolkit.displayln("CV parameters:");
        SFIToolkit.displayln("\t k = " + bestK);
        SFIToolkit.displayln("\t alpha = " + bestAlpha);
        SFIToolkit.displayln("\t mse_bar = " + bestMSEBar);
        SFIToolkit.displayln("\t depth = " + bestDepth);
        SFIToolkit.displayln("\t proportion used to estimate tree = " + 0.35);
        SFIToolkit.displayln("Best MSPE: " + bestMSPE);

        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * 0.35); // bestProportion);
        SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsToEstimateTreeStructure);
        Jama.Matrix treeX = new Jama.Matrix(numObsToEstimateTreeStructure, bartbootstrap.getX().getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsToEstimateTreeStructure, 1);
        Jama.Matrix honestX = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, bartbootstrap.getX().getColumnDimension());
        Jama.Matrix honestY = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, 1);
        TreeSet<Integer> treeSet = new TreeSet<>();
        int count = 0;

        Random honestRNG = new Random(787);

        while (count < numObsToEstimateTreeStructure) {
            if (count < numObsToEstimateTreeStructure) {
                int index = (int) Math.floor(honestRNG.nextDouble() * numObs);
                if (!treeSet.contains(index)) {
                    treeSet.add(index);
                    count++;
                }
            }
        }
        SFIToolkit.displayln("honestX row: " + honestX.getRowDimension() + " col: " + honestX.getColumnDimension());
        SFIToolkit.displayln("treeSet.size(): " + treeSet.size());

        int countTree = 0;
        int countHonest = 0;
        for (int i = 0; i < numObs; i++) {
            if (treeSet.contains(i)) {
                for (int j = 0; j < bartbootstrap.getX().getColumnDimension(); j++) {
                    treeX.set(countTree, j, bartbootstrap.getX().get(i, j));
                }
                treeY.set(countTree, 0, bartbootstrap.getY().get(i, 0));
                countTree++;

            } else {
                for (int j = 0; j < bartbootstrap.getX().getColumnDimension(); j++) {
                    honestX.set(countHonest, j, bartbootstrap.getX().get(i, j));
                }
                honestY.set(countHonest, 0, bartbootstrap.getY().get(i, 0));
                countHonest++;
            }
            // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
        }

//        System.out.println("treeX:");
//        pmUtility.prettyPrint(treeX);
//
//        System.out.println("honestX:");
//        pmUtility.prettyPrint(honestX);
        // this would be the place to introduce resampling of the treeX/Y and the honestX/Y to produce forest estimates
        ArrayList<TreeMoment> forest = new ArrayList<>();
        
        int numTrees = bartbootstrap.numberoftrees();
        // if this is too big, it doesn't run on Stata

        boolean verbose = false;
        boolean printTrees= false;

        Random rng = new Random();
        for (int fi = 0; fi < numTrees; fi++) {
            long seed = rng.nextLong();
            TreeMoment momentTree = new TreeMoment(null, bartbootstrap, resample(treeX, seed), resample(treeY, seed), bartbootstrap.getDiscreteVector(), verbose, bestAlpha, bestK, bestMSEBar, true, bestDepth);
            if (verbose) {
            	SFIToolkit.displayln("----------------------- mi: " + fi + " -----------------------");
            }
            momentTree.determineSplit();

            seed = rng.nextLong();
            Jama.Matrix thisHonestX = resample(honestX, seed);
            Jama.Matrix thisHonestY = resample(honestY, seed);

//            Jama.Matrix doublecheck = pmUtility.OLS(pmUtility.getColumn(thisHonestX, 0), thisHonestY, true);
//            pmUtility.prettyPrint(pmUtility.concatMatrix(thisHonestY, thisHonestX));
//            System.out.println("Manual: " + doublecheck.get(1, 0));
            for (int i = 0; i < thisHonestX.getRowDimension(); i++) {
                Jama.Matrix xi = thisHonestX.getMatrix(i, i, 0, thisHonestX.getColumnDimension() - 1);
                Jama.Matrix yi = thisHonestY.getMatrix(i, i, 0, thisHonestY.getColumnDimension() - 1);
                momentTree.sortXToCorrectLeafs(yi, xi);
            }
            momentTree.consolidateHonestData();
            momentTree.estimateHonestTree();
            if (verbose || printTrees) {
            	SFIToolkit.displayln("********* tree "+fi+" *********");
                momentTree.printTree();
            }
            forest.add(momentTree);
        }

        boolean writeOutput = true;

        if (writeOutput) {
            try {
                // BufferedWriter bartOut = new BufferedWriter(new FileWriter("bartTau.csv"));
                TreeSet<Double> uniqueTau = new TreeSet<>();
                TreeSet<Double> significantUniqueTau = new TreeSet<>();
                // Jama.Matrix beta_backto_Stata = new Jama.Matrix(allX.getRowDimension(),forest.size());
                int rcc;
                for (int i = 0; i < allXoriginal.getRowDimension(); i++) {
                    // for (int i = 0; i < 25; i++) {
                    Jama.Matrix xi = allXoriginal.getMatrix(i, i, 0, allXoriginal.getColumnDimension() - 1);              
                    xi.set(0, 0, 1); // set the treatment to one to see what the treatment effect is
                    Jama.Matrix tau = new Jama.Matrix(forest.size(), 1);
                    for (int mi = 0; mi < forest.size(); mi++) {
                        TreeMoment m = forest.get(mi);
                        Jama.Matrix beta = m.getEstimatedBeta(xi);
                        // pmUtility.prettyPrint(beta);
                        // tau.set(mi, 0, m.getPredictedY(xi));
                        // rcc  =  Data.storeNum(resv+1, i+1 ,beta.get(0, 0));
                        // Data.storeStr(resv+1, i ,String.valueOf(beta.get(0, 0)));
                        
                        if (verbose) {
                        	// SFIToolkit.displayln("mi: " + mi);
                        }
                        tau.set(mi, 0, beta.get(0, 0));
                    }               
                    uniqueTau.add(pmUtility.mean(tau, 0));
                    rcc  =  Data.storeNum(Data.getParsedVarCount()+2+resv, i+1 ,pmUtility.mean(tau, 0));
                    String stars = "";
                    NormalDistribution normal = new NormalDistribution();
                                       
//                 System.out.println(normal.inverse(0.05));
//                 System.out.println(normal.inverse(0.025));
//                 System.out.println(normal.inverse(0.005));
//                 System.exit(0);

                    boolean useZStat = false;
                    if (useZStat) {
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.05))) {
                            stars = "*";
                        }
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.025))) {
                            stars = "**";
                        }
                        if (Math.abs(pmUtility.mean(tau, 0) / pmUtility.standardDeviation(tau, 0)) > Math.abs(normal.inverse(0.005))) {
                            stars = "***";
                        }
                    } else {
                        if (pmUtility.percentile(tau, 0, 0.05) * pmUtility.percentile(tau, 0, 0.95) > 0) {
                            stars = "*";
                        }
                        if (pmUtility.percentile(tau, 0, 0.025) * pmUtility.percentile(tau, 0, 0.975) > 0) {
                            stars = "**";
                        }
                        if (pmUtility.percentile(tau, 0, 0.005) * pmUtility.percentile(tau, 0, 0.995) > 0) {
                            stars = "***";
                        }
                    }
                    if (!stars.equals("")) {
                        significantUniqueTau.add(pmUtility.mean(tau, 0));
                    }
                    System.out.format("%.3f \t %s \t %.3f \t [%.3f, %.3f] ", pmUtility.mean(tau, 0), stars, pmUtility.standardDeviation(tau, 0), pmUtility.percentile(tau, 0, 0.025), pmUtility.percentile(tau, 0, 0.975));
                    pmUtility.prettyPrint(xi);
                    // bartOut.write(pmUtility.mean(tau, 0) + "," + pmUtility.standardDeviation(tau, 0) + ",");
                    // for (int j = 0; j < xi.getColumnDimension(); j++) {
                    //    bartOut.write(xi.get(0, j) + ",");
                    // }
                    // bartOut.write("\n");
                }
                
                SFIToolkit.displayln("Estimating " + uniqueTau.size() + " unique treatment effects.");
                SFIToolkit.displayln("Estimating " + significantUniqueTau.size() + " statistically significant unique treatment effects.");
                // bartOut.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
         }
        

    }

    public static Jama.Matrix resample(Jama.Matrix x, long seed) {
        Random rng = new Random(seed);
        if (1 == 0) {
            return x.copy();
        }
        Jama.Matrix re = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension());
        for (int i = 0; i < re.getRowDimension(); i++) {
            int index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
            // SFIToolkit.displayln("INDEX " + index );
            for (int j = 0; j < re.getColumnDimension(); j++) {
                re.set(i, j, x.get(index, j));
            }
        }
        return re;
    }

}