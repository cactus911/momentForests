/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package examples.CardStata;

import Jama.Matrix;
import core.DataLens;
import core.HomogeneousSearchContainer;
import core.MomentForest;
import core.MomentSpecification;
import core.TreeMoment;
import core.TreeOptions;
import java.awt.BorderLayout;

import java.awt.GridLayout;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;
import java.util.StringJoiner;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.JTextAreaAutoscroll;
import utility.pmUtility;

import com.stata.sfi.Data;
import com.stata.sfi.SFIToolkit;
import com.stata.sfi.Macro;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
 */
public class CardMain {

    private ArrayList<Integer> homogeneousParameterList;
    private Jama.Matrix estimatedHomogeneousParameters;
    final private JTextArea jt;
    
    /**
     * @param args the command line arguments
     */
    
    public CardMain(JTextArea jt) {
        this.jt = jt;
    }
    
    public static void execute(CardSpecification spec, int seed) {
        // Redirect GUI output (JTextArea) to System.out for headless use
        JTextArea dummy = new JTextAreaAutoscroll() {
            @Override
            public void append(String str) {
                SFIToolkit.display(str);
            }
        };

        CardMain model = new CardMain(dummy);
        model.executeWithSpec(spec, seed);
    }

    public void executeWithSpec(CardSpecification spec, int seed) {
    	//Random rng = new Random(777);
        Random rng = new Random(seed);
        MomentSpecification mySpecification = spec;

        double bestMinImprovement = 4.0;
        int bestMinObservationsPerLeaf = 25;
        int bestMaxDepth = 0;

        double minOutOfSampleFit = 0;
        double minInSampleFit = 0;
        boolean first = true;
        
        /**
         * Make sure that cross-validation is run on completely unrestricted
         * model; set all parameters to heterogeneous
         */
        mySpecification.resetHomogeneityIndex();

        int numberTreesInForest = spec.getNumTrees();
        double proportionObservationsToEstimateTreeStructure = spec.getProportionObservationsToEstimateTreeStructure();
        boolean detectHomogeneity = spec.doDetectHomogeneity();
        // SFIToolkit.displayln("numTrees: " + numberTreesInForest);

        /*
         * Initialize the moment forest
         */
        boolean verbose = false;
        long rngBaseSeedMomentForest = rng.nextLong();
        long rngBaseSeedOutOfSample = rng.nextLong();
        
        boolean runCV = spec.doCrossValidation();
        if (runCV) {
            SFIToolkit.displayln("****************************");
            SFIToolkit.displayln("* Running Cross-Validation *");
            SFIToolkit.displayln("****************************");

            ArrayList<computeFitStatistics> cvList = new ArrayList<>();
            List<Integer> gridMinLeaf = spec.getCVGridMinLeaf();
            List<Double> gridMinImprovement = spec.getCVGridMinImprovement();
            List<Integer> gridMaxDepth = spec.getCVGridMaxDepth();
            
            int[] stratIndices = mySpecification.getStratificationIndex();
            if (stratIndices != null && stratIndices.length > 0) {
                StringJoiner sj = new StringJoiner(", ");
                for (int idx : stratIndices) {
                    sj.add(mySpecification.getVariableName(idx));
                }
                SFIToolkit.displayln(
                    "Splitting sample using stratified random sampling on variables: " + sj.toString()
                );
            } else {
            	SFIToolkit.displayln("Splitting sample using simple random sampling.");
            }
            
            /*
            for (int minObservationsPerLeaf = 25; minObservationsPerLeaf <= 200; minObservationsPerLeaf *= 2) {
                for (double minImprovement = 0.1; minImprovement <= 2.0; minImprovement *= 2) {
                    for (int maxDepth = 7; maxDepth >= 1; maxDepth--) {
                    	cvList.add(new computeFitStatistics(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample, false));
                    }
                }
            }
            */
       
            SFIToolkit.displayln("Grid for minLeaf: " + gridMinLeaf.toString());
            SFIToolkit.displayln("Grid for minImprovement: " + gridMinImprovement.toString());
            SFIToolkit.displayln("Grid for maxDepth: " + gridMaxDepth.toString());
            SFIToolkit.displayln("Proportion of observations to estimate tree structure: " + proportionObservationsToEstimateTreeStructure);
                             
            for (int minObservationsPerLeaf : gridMinLeaf) {
                for (double minImprovement : gridMinImprovement) {
                    for (int maxDepth : gridMaxDepth) {
                    	
                    	
                    	//SFIToolkit.displayln("Parameters before detection (minLeaf=" + minObservationsPerLeaf + ", minImprovement=" + minImprovement + ", maxDepth=" + maxDepth + "): " + java.util.Arrays.toString(mySpecification.getHomogeneousIndex()));                	
                    	mySpecification.resetHomogeneityIndex();
                        if (detectHomogeneity) {
                            executeHomogeneousParameterClassificationAndSearch(mySpecification, numberTreesInForest, minImprovement, minObservationsPerLeaf, maxDepth, rngBaseSeedMomentForest, rngBaseSeedOutOfSample, false, false);
                        }                      
                        //SFIToolkit.displayln("Parameters flags after detection (minLeaf=" + minObservationsPerLeaf + ", minImprovement=" + minImprovement + ", maxDepth=" + maxDepth + "): " + java.util.Arrays.toString(mySpecification.getHomogeneousIndex()));
                        computeFitStatistics s = new computeFitStatistics(mySpecification, numberTreesInForest, proportionObservationsToEstimateTreeStructure, rngBaseSeedMomentForest, verbose, minObservationsPerLeaf, minImprovement, maxDepth, rngBaseSeedOutOfSample, false);
                        s.computeOutOfSampleMSE();
                        cvList.add(s);
                    }
                }
            }
            
            /*
            cvList.parallelStream().forEach(s -> {
                s.computeOutOfSampleMSE();
            });
			*/
                 
            for (computeFitStatistics s : cvList) {
                double combinationMSE = s.getMSE();
                String star = "";
                String starIn = "";
                if (combinationMSE <= minOutOfSampleFit || first) {
                    minOutOfSampleFit = combinationMSE;
                    bestMinImprovement = s.getMinImprovement();
                    bestMinObservationsPerLeaf = s.getMinObservationsPerLeaf();
                    bestMaxDepth = s.getMaxTreeDepth();
                    star = "(*)";
                }
                //SFIToolkit.display("best: " + String.format("%.6f", minInSampleFit) + " this: " + String.format("%.6f", s.getInSampleFit()) + " " + (minInSampleFit > s.getInSampleFit()));
                if (s.getInSampleFit() < minInSampleFit) {
                	if (verbose) {
                		SFIToolkit.displayln("detected lower");
                	}
                    minInSampleFit = s.getInSampleFit();
                    starIn = "(**)";
                }

                if (first) {
                    first = false;
                }
                if (verbose) {
                	SFIToolkit.displayln(starIn);
                	SFIToolkit.displayln("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + String.format("%.10f", combinationMSE) + " " + star);
                	//jt.append("minMSE: " + s.getMinImprovement() + " minObs: " + s.getMinObservationsPerLeaf() + " maxDepth: " + s.getMaxTreeDepth() + " Out-of-sample MSE: " + combinationMSE + " " + star + " In-Sample Fit: " + s.getInSampleFit() + " " + starIn + "\n");
                }
            }
            SFIToolkit.displayln("Lowest MSE: " + String.format("%.10f", minOutOfSampleFit) + " at min_N = " + bestMinObservationsPerLeaf + " min_improvement = " + bestMinImprovement + " maxDepth: " + bestMaxDepth);
            //jt.append("Lowest MSE: " + minOutOfSampleFit + " at min_N = " + bestMinObservationsPerLeaf + " min_MSE = " + bestMinImprovement + " maxDepth: " + bestMaxDepth + "\n");
            //jt.append("Best in-sample fit: " + minInSampleFit + "\n");
        } else {
            bestMinObservationsPerLeaf = 25;
            bestMinImprovement = 0.1;
            bestMaxDepth = 5;
        }
                           
        if (detectHomogeneity && bestMaxDepth > 0) {
            executeHomogeneousParameterClassificationAndSearch(mySpecification, numberTreesInForest, bestMinImprovement, bestMinObservationsPerLeaf, bestMaxDepth, rngBaseSeedMomentForest, rngBaseSeedOutOfSample, false, true);
        }
        
        /**
         * Compute out-of-sample measures of fit (against Y, and true beta)
         */      
        SFIToolkit.displayln("*******************************************");
        SFIToolkit.displayln("* Computing out-of-sample measures of fit *");
        SFIToolkit.displayln("*******************************************");
        
        computeFitStatistics fitStats = new computeFitStatistics(mySpecification, numberTreesInForest, proportionObservationsToEstimateTreeStructure, rngBaseSeedMomentForest, verbose, bestMinObservationsPerLeaf,
                bestMinImprovement, bestMaxDepth, rngBaseSeedOutOfSample, true);
        fitStats.computeOutOfSampleMSE();
        double outOfSampleFit = fitStats.getMSE();

        SFIToolkit.displayln("Best minimum observations in leaf: " + bestMinObservationsPerLeaf);
        SFIToolkit.displayln("Best minimum improvement: " + bestMinImprovement);
        SFIToolkit.displayln("Best maximum depth: " + bestMaxDepth);

        SFIToolkit.displayln("Out of sample SSE: " + String.format("%.10f", outOfSampleFit));
        Macro.setLocal("mf_oos_fit", String.valueOf(outOfSampleFit));
        
        SFIToolkit.displayln("Number of trees in forest: " + numberTreesInForest);
        
        double[] countVariableSplitsInForest = fitStats.getSplitVariables();
        //SFIToolkit.displayln("Number of split variables: " + countVariableSplitsInForest.length);

        SFIToolkit.displayln("Forest split on the following variables:");
        for (int i = 0; i < countVariableSplitsInForest.length; i++) {
            if (countVariableSplitsInForest[i] > 0) {
                String varName = mySpecification.getVariableName(i);
                double pct = 100.0 * countVariableSplitsInForest[i] / numberTreesInForest;
                String line = String.format("%20s [%.2f%%]", varName, pct);
                SFIToolkit.displayln(line);
            }
        }
        
        fitStats.outputInSampleFits();
    }
      
    private void executeHomogeneousParameterClassificationAndSearch(MomentSpecification mySpecification, int numberTreesInForest, double minImprovement, int minObservationsPerLeaf, int maxDepth, long rngBaseSeedMomentForest, long rngBaseSeedOutOfSample, boolean verbose, boolean verboselast) {        	
        
    	if (verboselast) {
        	SFIToolkit.displayln("***************************");
            SFIToolkit.displayln("* Testing for Homogeneity *");
            SFIToolkit.displayln("***************************");
        }

        double minProportionInEachLeaf = 0.01;

        // Create forest over current specification
        DataLens forestLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());
        TreeOptions cvOptions = new TreeOptions(minProportionInEachLeaf, minObservationsPerLeaf, minImprovement, maxDepth, true); 
        MomentForest myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, forestLens, verbose, new TreeOptions());
        
        myForest.setTreeOptions(cvOptions);
        myForest.growForest();
            
        if (verbose) {
            TreeMoment loblolly = myForest.getTree(0);
            SFIToolkit.displayln("************************");
            SFIToolkit.displayln("* Printing first tree  *");
            SFIToolkit.displayln("************************");
            loblolly.printTree();
        }
           
        myForest.testHomogeneity(verboselast);
        
        // Collect homogeneity votes and starting values
        ArrayList<Integer> hpl = new ArrayList<>();
        ArrayList<Double> hplStartingValues = new ArrayList<>();
        boolean[] voteIndexHomogeneity = myForest.getHomogeneityVotes(jt, verboselast);
        
        double[] startingValues = myForest.getHomogeneityStartingValues();
        for (int i = 0; i < voteIndexHomogeneity.length; i++) {
            if (voteIndexHomogeneity[i]) {
            	if (verboselast) {
            		SFIToolkit.displayln("Adding index " + i + " to homogeneous list with starting value: " + startingValues[i]);
            	}
                hpl.add(i);
                hplStartingValues.add(startingValues[i]);
            }
        }

        // Sort and set indices
        HomogeneousParameterSorter sorter = new HomogeneousParameterSorter();
        sorter.sort(hpl, hplStartingValues);
        
        mySpecification.resetHomogeneityIndex();
        for (int i = 0; i < hpl.size(); i++) {
            mySpecification.setHomogeneousIndex(hpl.get(i));
            mySpecification.setHomogeneousParameter(hpl.get(i), hplStartingValues.get(i));
        }
        setHomogeneousParameterList(hpl);
        boolean testPostClassificationConvergence = !true;
        if (testPostClassificationConvergence) {
            hpl.clear();
            hpl.add(0);
            mySpecification.resetHomogeneityIndex();
            mySpecification.setHomogeneousIndex(0);
            mySpecification.setHomogeneousParameter(0, -1.0);
        }
        
        
        /*
         * Estimate values of those homogeneous parameters
         */
        if (!hpl.isEmpty()) {
            if (hpl.size() == mySpecification.getNumParams()) {
                // all homogenous, don't need to optimize, just set to stump and let it run
            	maxDepth = 0;
            } else {
                SFIToolkit.displayln("Initializing search container");
                HomogeneousSearchContainer con = new HomogeneousSearchContainer(mySpecification, numberTreesInForest, verbose, minImprovement, minObservationsPerLeaf, maxDepth, getHomogeneousParameterList(), rngBaseSeedMomentForest, rngBaseSeedOutOfSample);
                SFIToolkit.displayln("Calling execute search");
                con.executeSearch();
                SFIToolkit.displayln("Post search");
                
                Jama.Matrix homogeneousParameters = con.getEstimatedHomogeneousParameters();
                SFIToolkit.display("Post-HomogeneousSearchContainer Estimated homogeneous parameters: ");
                pmUtility.prettyPrintVector(homogeneousParameters);

                int K = mySpecification.getHomogeneousParameterVector().getRowDimension();
                SFIToolkit.displayln("Post-HomogeneousSearchContainer Length of homogeneous parameter vector: " + K);
                Jama.Matrix expandedHomogeneousParameterVector = new Jama.Matrix(K, 1);
                int counter = 0;
                for (int k = 0; k < K; k++) {
                    if (mySpecification.getHomogeneousIndex()[k]) {
                        expandedHomogeneousParameterVector.set(k, 0, homogeneousParameters.get(counter, 0));
                        mySpecification.setHomogeneousParameter(k, homogeneousParameters.get(counter, 0));
                        counter++;
                    }
                }
                SFIToolkit.displayln("Specification homogeneous parameter vector: ");
                pmUtility.prettyPrintVector(mySpecification.getHomogeneousParameterVector());
                // System.exit(0);
                setEstimatedHomogeneousParameters(expandedHomogeneousParameterVector);
            }
        }
    }
    
    private class computeFitStatistics {

        MomentSpecification mySpecification;
        int numberTreesInForest;
        double proportionObservationsToEstimateTreeStructure;
        long rngBaseSeedMomentForest;
        boolean verbose;
        int minObservationsPerLeaf;
        double minImprovement;
        int maxTreeDepth;
        long rngBaseSeedOutOfSample;
        boolean generatePlots;

        double MSE;
        double inSampleFit;

        MomentForest myForest;

        public computeFitStatistics(MomentSpecification mySpecification, int numberTreesInForest, double proportionObservationsToEstimateTreeStructure, long rngBaseSeedMomentForest, boolean verbose, int minObservationsPerLeaf, double minImprovement, int maxTreeDepth, long rngBaseSeedOutOfSample, boolean generatePlots) {
            this.mySpecification = mySpecification;
            this.numberTreesInForest = numberTreesInForest;
            this.proportionObservationsToEstimateTreeStructure = proportionObservationsToEstimateTreeStructure;
            this.rngBaseSeedMomentForest = rngBaseSeedMomentForest;
            this.verbose = verbose;
            this.minObservationsPerLeaf = minObservationsPerLeaf;
            this.minImprovement = minImprovement;
            this.maxTreeDepth = maxTreeDepth;
            this.rngBaseSeedOutOfSample = rngBaseSeedOutOfSample;
            this.generatePlots = generatePlots;
        }

        public double getMSE() {
            return MSE;
        }

        public int getMaxTreeDepth() {
            return maxTreeDepth;
        }

        public double getMinImprovement() {
            return minImprovement;
        }

        public int getMinObservationsPerLeaf() {
            return minObservationsPerLeaf;
        }

        public double[] getSplitVariables() {
            return myForest.getNumberTimesTreesInForestSplitOnAGivenVariableIndex();
        }

        public void computeOutOfSampleMSE() {
            // SFIToolkit.displayln("\nComputing OOS In Parameter Space\n");
            // SFIToolkit.displayln("Homogeneous parameter length in spec: "+mySpecification.getHomogeneousIndex().length);
            DataLens overallLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());

            DataLens[] split;
            if (mySpecification.getStratificationIndex() != null) {
            	//SFIToolkit.displayln("Splitting sample using stratified random sampling on column index: " + mySpecification.getStratificationIndex());
                split = overallLens.randomlySplitSampleByStrata(0.9, rngBaseSeedMomentForest);
            } else {
            	//SFIToolkit.displayln("Splitting sample using simple random sampling.");
            	split = overallLens.randomlySplitSample(0.9, rngBaseSeedMomentForest);
            }
            
            DataLens estimatingLens = split[0];
            DataLens oosDataLens = split[1];
            
            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, estimatingLens, verbose, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();
            if (verbose) {
                SFIToolkit.displayln("First tree in forest estimated as:");
                myForest.getTree(0).printTree();
            }
            /**
             * Test vectors for assessment
             */
            // DataLens oosDataLens = mySpecification.getOutOfSampleXYZ(1500, rngBaseSeedOutOfSample); // this should eventually be modified to come out of the data itself (or generalized in some way)
            Jama.Matrix testZ = oosDataLens.getZ();
            Jama.Matrix testX = oosDataLens.getX();
            Jama.Matrix testY = oosDataLens.getY();

            double outOfSampleFit = 0;
            for (int i = 0; i < testZ.getRowDimension(); i++) {
                Jama.Matrix zi = testZ.getMatrix(i, i, 0, testZ.getColumnDimension() - 1);
                Jama.Matrix xi = testX.getMatrix(i, i, 0, testX.getColumnDimension() - 1);
                double yi = testY.get(i, 0);

                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
                outOfSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
            }

            inSampleFit = 0;
            for (int i = 0; i < mySpecification.getZ().getRowDimension(); i++) {
                Jama.Matrix zi = mySpecification.getZ().getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);
                Jama.Matrix xi = mySpecification.getX().getMatrix(i, i, 0, mySpecification.getX().getColumnDimension() - 1);
                double yi = mySpecification.getY().get(i, 0);
                // have to reconstruct a composite beta from homogeneous and heterogeneous parameters
                Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);

                inSampleFit += mySpecification.getGoodnessOfFit(yi, xi, compositeEstimatedBeta);
            }
            inSampleFit /= mySpecification.getZ().getRowDimension();

            MSE = outOfSampleFit / testZ.getRowDimension(); // mse
     
        }

        public void outputInSampleFits() {
            boolean verboseInSample = false;
            DataLens overallLens = new DataLens(mySpecification.getX(), mySpecification.getY(), mySpecification.getZ(), null, mySpecification.getStratificationIndex());

            myForest = new MomentForest(mySpecification, numberTreesInForest, rngBaseSeedMomentForest, overallLens, verboseInSample, new TreeOptions());
            TreeOptions cvOptions = new TreeOptions(0.01, minObservationsPerLeaf, minImprovement, maxTreeDepth, false); // k = 1
            myForest.setTreeOptions(cvOptions);
            /**
             * Grow the moment forest
             */
            myForest.growForest();
            SFIToolkit.displayln("First tree in forest estimated as:");
            myForest.getTree(0).printTree();
            
            try {
            	if(mySpecification.getBetaPrefixes() != "") {
            		
            		// Create variables in Stata dataset
            		for(int j=0; j < mySpecification.getX().getColumnDimension(); j++) {
            			String betaVarName = mySpecification.getBetaPrefixes() + j;
            			if (Data.getVarIndex(betaVarName) != 0) {
            				throw new IllegalStateException("Variable " + betaVarName + " already exists.");   
            			} else {
            				Data.addVarDouble(betaVarName);
            			}
            		}
            		
            		// Store beta values into Stata
            		for (int i = 0; i < Data.getObsTotal(); i++) {
            			Jama.Matrix zi = mySpecification.getZ().getMatrix(i, i, 0, mySpecification.getZ().getColumnDimension() - 1);
            			Jama.Matrix compositeEstimatedBeta = myForest.getEstimatedParameterForest(zi);
            			
                        for(int j=0; j < mySpecification.getX().getColumnDimension(); j++) {
                        	String betaVar = mySpecification.getBetaPrefixes() + j;
                        	Data.storeNum(Data.getVarIndex(betaVar), i+1, compositeEstimatedBeta.get(j,0));
                        }
            		}
            	}
            	
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(0);
            }

        }

        private double getInSampleFit() {
            return inSampleFit;
        }
    }

    private void setHomogeneousParameterList(ArrayList<Integer> homogeneousParameterList) {
        this.homogeneousParameterList = homogeneousParameterList;
    }

    public ArrayList<Integer> getHomogeneousParameterList() {
        return homogeneousParameterList;
    }

    /**
     * @return the estimateHomogeneousParameter
     */
    public Jama.Matrix getEstimatedHomogeneousParameters() {
        return estimatedHomogeneousParameters;
    }

    public void setEstimatedHomogeneousParameters(Matrix estimatedHomogeneousParameters) {
        this.estimatedHomogeneousParameters = estimatedHomogeneousParameters;
    }

}
