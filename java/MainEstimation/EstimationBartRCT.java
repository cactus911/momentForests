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
package MainEstimation;

import com.stata.sfi.*;
import JSci.maths.statistics.NormalDistribution;
import MainEstimation.TreeMoment;
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
public class EstimationBartRCT {
	
	 public static int Momentforests(String[] args) { 
		 
     // public static void main(String[] args) {     

		int resv = Data.getParsedVarCount(); // Data.getVarIndex(args[0]);  
		
		// SFIToolkit.displayln("Boostrap Count: " + resv);
                int rc = 0;
		EstimationBartRCT go = new EstimationBartRCT(resv);
		
		return(rc);
     }
    
	
	// public static int loadData(String[] args) { 
       @SuppressWarnings("deprecation")
	public EstimationBartRCT(int resv) {
   		MomentSpecificationBartRCT bart = new MomentSpecificationBartRCT();
	    bart.loadData();		 
       
        int bestK = 0;
        double bestAlpha = 0;
        double bestMSEBar = 0;
        int bestDepth = 0;
        boolean first = true;
        double bestMSPE = 0;
        double bestProportion = 0;
        
	    int rc_cv = 0 ;

        int numObs = bart.getX().getRowDimension();
        Jama.Matrix allX = bart.getX().copy();
 
        // try implementing a true OOS prediction error here?
        for (double proportion = 0.3; proportion <= 0.3; proportion += 0.05) {
            int numObsEstimateLeafValues = (int) Math.floor(numObs * proportion);
            int numPredictObs = (int) Math.floor(numObs * 0.1); // it was originally 10;
            int numObsGrowTreeStructure = numObs - numObsEstimateLeafValues - numPredictObs; // halfObs / 2;
            SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsEstimateLeafValues + " predictObs: " + numPredictObs);

            Jama.Matrix treeX = new Jama.Matrix(numObsGrowTreeStructure, bart.getX().getColumnDimension());
            Jama.Matrix treeY = new Jama.Matrix(numObsGrowTreeStructure, 1);
            Jama.Matrix honestX = new Jama.Matrix(numObsEstimateLeafValues, bart.getX().getColumnDimension());
            Jama.Matrix honestY = new Jama.Matrix(numObsEstimateLeafValues, 1);
            Jama.Matrix predictX = new Jama.Matrix(numPredictObs, bart.getX().getColumnDimension());
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
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        treeX.set(countTree, j, bart.getX().get(i, j));
                    }
                    treeY.set(countTree, 0, bart.getY().get(i, 0));  
                    countTree++;
                } else if (predictSet.contains(i)) {
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        predictX.set(countPredict, j, bart.getX().get(i, j));
                    }
                    predictY.set(countPredict, 0, bart.getY().get(i, 0));
                    countPredict++;
                } else {
                    for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                        honestX.set(countHonest, j, bart.getX().get(i, j));
                    }
                    honestY.set(countHonest, 0, bart.getY().get(i, 0));
                    countHonest++;
                }
                // System.out.println(i+" "+treeSet.contains(i)+" "+countTree+" "+countHonest);
            }

            double minProportionEachPartition = 0.001;
            int minCountEachPartition = 10;
            double improvementThreshold = 0.01;
            int maxDepth = 100;
            boolean verbose = false;
            
            Jama.Matrix CVParameters = bart.cvparameters();
            
            int minCountEachPartition_start = (int) CVParameters.get(0, 0);
            int minCountEachPartition_jumpsize = (int) CVParameters.get(0, 1);
            int minCountEachPartition_end = (int) CVParameters.get(0, 2);
            
            
            double improvementThreshold_start = CVParameters.get(0, 3); 
            double improvementThreshold_jumpsize = CVParameters.get(0, 4);
            double improvementThreshold_end = CVParameters.get(0, 5);

            
            int dim1 = (minCountEachPartition_end-minCountEachPartition_start)/minCountEachPartition_jumpsize;
            double dim2 = (improvementThreshold_end-improvementThreshold_start)/improvementThreshold_jumpsize;
            int dim22 = (int) dim2;
            int Dim = (int) (dim1+1)*(dim22+1);
            
            Jama.Matrix CV_parameters = new Jama.Matrix(Dim, 3);
            
            
            // add randomization layer here
            Random rng_cv1 = new Random();
            long seed_cv1 = rng_cv1.nextLong();
            
            Jama.Matrix mink;
            mink = new Jama.Matrix(dim1, 1);
            int mink_value = minCountEachPartition_start;
            for (int i = 0; i < dim1; i++) {
            	mink.set(i, 0, mink_value);
            	mink_value = mink_value + minCountEachPartition_jumpsize;
            }
            
            Jama.Matrix mink_random = resample(mink, seed_cv1);
            
            Random rng_cv2 = new Random();
            long seed_cv2 = rng_cv1.nextLong();
            
            Jama.Matrix barmse;
            
            barmse = new Jama.Matrix(dim22, 1);
            double barmse_value = improvementThreshold_start;
            for (int i = 0; i < dim22; i++) {
            	barmse.set(i, 0, barmse_value);
            	barmse_value = barmse_value + improvementThreshold_jumpsize;
            }
            
            Jama.Matrix barmse_random = resample(barmse, seed_cv2);
           

            int CV_index = 1;
            // for (maxDepth = 1; maxDepth <= 20; maxDepth++) {
            for (int ii = 0; ii < dim1; ii++) {         
                minCountEachPartition = (int) mink_random.get(ii,0);
            // for (minCountEachPartition = minCountEachPartition_end /*40*/; minCountEachPartition >= minCountEachPartition_start; minCountEachPartition -= minCountEachPartition_jumpsize) {
                for (int j = 0; j < dim22; j++) {
                	improvementThreshold = barmse_random.get(j,0);
                // for (improvementThreshold = improvementThreshold_end /*0.07514298*/; improvementThreshold >= improvementThreshold_start; improvementThreshold -= improvementThreshold_jumpsize) {
                	SFIToolkit.displayln("MSE_bar: " + improvementThreshold + " k: " + minCountEachPartition);
                    TreeMoment momentTree = new TreeMoment(null, bart, treeX, treeY, bart.getDiscreteVector(), verbose, minProportionEachPartition, minCountEachPartition, improvementThreshold, true, maxDepth);
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
                    int CV_Index = CV_index - 1;
                    SFIToolkit.displayln("CV_index: " + CV_Index + " " );
                    SFIToolkit.displayln("maxDepth: " + maxDepth + " " + " k: " + minCountEachPartition + " " + " alpha: " + minProportionEachPartition + " " + " mse_bar: " + improvementThreshold + " ");
                    SFIToolkit.displayln("MSE: " + MSE + " MSPE: " + MSPE + " nulls: " + nullCounter + " ");
                    // System.out.println("MSPE: " + MSPE + " nulls: " + nullCounterMSPE);
                    if (MSPE < bestMSPE || first) {
                    	SFIToolkit.displayln(" * ");
                        // bestK = minCountEachPartition;
                        bestAlpha = minProportionEachPartition;
                        // bestMSEBar = improvementThreshold;
                        bestDepth = maxDepth;
                        // bestMSPE = MSPE;
                        bestProportion = proportion;
                        first = false;
                    } else {
                    	SFIToolkit.displayln(" ");
                    }
                    
                    CV_parameters.set(CV_index-1, 0, minCountEachPartition);
                    CV_parameters.set(CV_index-1, 1, improvementThreshold);
                    CV_parameters.set(CV_index-1, 2, MSPE);
                    
                    CV_index ++;                

                }
            }
            
    		for (int i = 0; i < Dim; i++) {
    			System.out.println(" Min Index : " + i + " Min SSPE : " + CV_parameters.get(i,2)  + " Min k : " + CV_parameters.get(i,0) + " Min MSE_bar : " + CV_parameters.get(i,1) ); 
    				// minAt = CV_parameters.get(i,2) > CV_parameters.get(minAt,2) ? i : minAt;
    		}
    		
            // Find out the index of those minimizing CV parameters 
            int minAt = 0;
    		for (int i = 1; i < Dim; i++) {
    			 if ( CV_parameters.get(i,2) < CV_parameters.get(minAt,2) && CV_parameters.get(i,2) > 0) minAt = i; 
    				// minAt = CV_parameters.get(i,2) > CV_parameters.get(minAt,2) ? i : minAt;
    		}
    		System.out.println(" Found Min Index : " + minAt + " Min SSPE : " + CV_parameters.get(minAt,2)  + " Min k : " + CV_parameters.get(minAt,0) + " Min MSE_bar : " + CV_parameters.get(minAt,1) );		  
           
    		bestMSPE = CV_parameters.get(minAt,1);
    		bestMSEBar =  CV_parameters.get(minAt,1);  //  
    		bestK = (int) CV_parameters.get(minAt,0);
    		
    		
            rc_cv = Data.storeNum(1, numObs+3, bestK);
            rc_cv = Data.storeNum(1, numObs+4, bestMSEBar);
            
        }
        SFIToolkit.displayln("CV parameters:");
        SFIToolkit.displayln("\t k = " + bestK);
        SFIToolkit.displayln("\t alpha = " + bestAlpha);
        SFIToolkit.displayln("\t mse_bar = " + bestMSEBar);
        SFIToolkit.displayln("\t depth = " + bestDepth);
        SFIToolkit.displayln("\t proportion used to estimate tree = " + 0.35); // bestProportion);
        SFIToolkit.displayln("Best MSPE: " + bestMSPE);
        
        // Let's put these CV parameters separately so that the same CV parameters would be used for bootstrapping


        int numObsToEstimateTreeStructure = (int) Math.floor(numObs * 0.35);
        SFIToolkit.displayln("numObs: " + numObs + " halfObs: " + numObsToEstimateTreeStructure);
        Jama.Matrix treeX = new Jama.Matrix(numObsToEstimateTreeStructure, bart.getX().getColumnDimension());
        Jama.Matrix treeY = new Jama.Matrix(numObsToEstimateTreeStructure, 1);
        Jama.Matrix honestX = new Jama.Matrix(numObs - numObsToEstimateTreeStructure, bart.getX().getColumnDimension());
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
                for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                    treeX.set(countTree, j, bart.getX().get(i, j));
                }
                treeY.set(countTree, 0, bart.getY().get(i, 0));
                countTree++;

            } else {
                for (int j = 0; j < bart.getX().getColumnDimension(); j++) {
                    honestX.set(countHonest, j, bart.getX().get(i, j));
                }
                honestY.set(countHonest, 0, bart.getY().get(i, 0));
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
        
        int numTrees = bart.numberoftrees();
        // int numTrees = 1; // 500; if this is too big, it doesn't run on Stata

        boolean verbose = true;
        boolean printTrees= false;

        Random rng = new Random();
        for (int fi = 0; fi < numTrees; fi++) {
            long seed = rng.nextLong();
            TreeMoment momentTree = new TreeMoment(null, bart, resample(treeX, seed), resample(treeY, seed), bart.getDiscreteVector(), verbose, bestAlpha, bestK, bestMSEBar, true, bestDepth);
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
                for (int i = 0; i < allX.getRowDimension(); i++) {
                    // for (int i = 0; i < 25; i++) {
                    Jama.Matrix xi = allX.getMatrix(i, i, 0, allX.getColumnDimension() - 1);              
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
                        // beta_backto_Stata.set(i, mi, beta.get(0, 0));
                    }               
                    uniqueTau.add(pmUtility.mean(tau, 0));
                    // SFIToolkit.displayln("Is this real???" + pmUtility.mean(tau, 0));
                    rcc  =  Data.storeNum(resv+2, i+1 ,pmUtility.mean(tau, 0));
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
                // SFIToolkit.displayln("beta Row Dim" + beta_backto_Stata.getRowDimension());  
                // SFIToolkit.displayln("beta Col Dim" + beta_backto_Stata.getColumnDimension()); 
                
              //Create variables in order to copy results back to Stata
                int rc;
                // rc = Data.addVarLong("beta_estimated");
                // if (rc!=0) return(rc); 
                
                
                //private static void processData(long num_rows, beta_backto_Stata){
                    // int                     rc;
                	// SFIToolkit.displayln("beta HERE" + beta_backto_Stata.get(0,0)); 
                	// SFIToolkit.displayln("resv: " + resv);
                	SFIToolkit.displayln("getParsedVarCount: " + Data.getParsedVarCount() );
                    // if (rc!=0) return(rc);

                    // return(rc);
           // }
                              
                
                SFIToolkit.displayln("Estimating " + uniqueTau.size() + " unique treatment effects.");
                SFIToolkit.displayln("Estimating " + significantUniqueTau.size() + " statistically significant unique treatment effects.");
                // bartOut.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
         }
        
    // return(rc);
    
    }

    public static Jama.Matrix resample(Jama.Matrix x, long seed) {
        Random rng = new Random(seed);
        if (1 == 0) {
            return x.copy();
        }
        Jama.Matrix re = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension());
        for (int i = 0; i < re.getRowDimension(); i++) {
            int index = (int) Math.floor(re.getRowDimension() * rng.nextDouble());
            for (int j = 0; j < re.getColumnDimension(); j++) {
                re.set(i, j, x.get(index, j));
            }
        }
        return re;
    }
}
