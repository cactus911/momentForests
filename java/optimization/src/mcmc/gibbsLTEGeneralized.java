/*
 * gibbsLTE.java
 *
 * Created on May 10, 2004, 2:39 PM
 */
package mcmc;

import JSci.maths.statistics.*;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import javax.swing.*;
import javax.swing.table.TableModel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan
 */
public class gibbsLTEGeneralized {

    public boolean gui = true;
    boolean TUNE_VARIANCE = false;
    boolean GENERATE_SET_SAMPLES = false;
    boolean ACCEPT_ONLY_IMPROVEMENTS = false;  // this is if you want to only accept improvements
    boolean DRAW_ALL_SIMULTANEOUSLY = false;
    int STOPPING_CRITERION = 250;
    boolean USE_ADAPTIVE_STOPPING_RULE = false;
    boolean quiet = true; // do not report any information to the command line or use a gui
    Jama.Matrix chain;
    double acceptance_ratio = 0.234;
    double minimumVariance = 0.1;
    double[] variance;
    double[] savedVar;
    java.text.NumberFormat nf = java.text.NumberFormat.getInstance();
    boolean loadedVar = false;
    /*
     * GUI elements; do not initialize unless using the gui to prevent problems
     * when LTE is called many times.
     */
    javax.swing.JLabel fval;
    javax.swing.JLabel[] labels;
    javax.swing.JProgressBar bar;
    javax.swing.JFrame f;
    JCheckBox showParameters;
    JCheckBox showValue;
    JCheckBox keepGoing;
    JCheckBox logYAxis;
    JButton executeNowButton;
    JButton exitNow;
    JButton shrinkVarianceButton;
    JTable statusTable;
    JLabel movingAverageLabel;    // add a graph for best value and current value
    JFreeChart chart;
    String[] varLabels;
    private int burnProgress = 0;
    boolean haveLabels = false;
    private double[] lowestPoint;
    private double[] point;
    double lowestValue;
    XYSeriesCollection xyc = new XYSeriesCollection();
    XYSeries xyBest = new XYSeries("Best");
    XYSeries xyCurrent = new XYSeries("Current");
    static int NUM_INSTANCES = 0;
    boolean continueExecution = true;
    private double TARGET_THRESHOLD = 0;
    private boolean USE_TARGET_THRESHOLD = false;
    private boolean[] ignoreList = null;

    String logFilename = null;

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String varFile, String[] varLabels) {
        this.varLabels = varLabels;
        haveLabels = true;
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream("data/" + varFile + ".dat"));
            variance = new double[guess.length];
            savedVar = new double[guess.length];
            for (int i = 0; i < guess.length; i++) {
                variance[i] = in.readDouble();
                savedVar[i] = variance[i];
            }
            in.close();
            // System.out.println("Read in variance vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(variance, 1)));
            loadedVar = true;
        } catch (Exception e) {
            // e.printStackTrace();
            // System.out.println("Creating new variance storage file: "+varFile+".dat");
        }
        execute(function, draws, burnPeriod, guess);
        try {
            ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream("data/" + varFile + ".dat"));
            for (int i = 0; i < guess.length; i++) {
                out.writeDouble(variance[i]);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String varFile) {
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream("data/" + varFile + ".dat"));
            variance = new double[guess.length];
            savedVar = new double[guess.length];
            for (int i = 0; i < guess.length; i++) {
                variance[i] = in.readDouble();
                savedVar[i] = variance[i];
            }
            in.close();
            // System.out.println("Read in variance vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(variance, 1)));
            loadedVar = true;
        } catch (Exception e) {
            // e.printStackTrace();
            System.out.println("Creating new variance storage file: " + varFile + ".dat");
        }
        execute(function, draws, burnPeriod, guess);
        try {
            ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream("data/" + varFile + ".dat"));
            for (int i = 0; i < guess.length; i++) {
                out.writeDouble(variance[i]);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess) {
        execute(function, draws, burnPeriod, guess);
    }

    /**
     * Creates a new instance of gibbsLTE
     */
    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels) {
        this.varLabels = varLabels;
        haveLabels = true;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, boolean useGUI) {
        gui = useGUI;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, boolean useGUI, boolean[] ignoreListPassed) {
        gui = useGUI;
        this.ignoreList = ignoreListPassed;

        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, boolean useGUI, boolean useAdaptive) {
        gui = useGUI;
        USE_ADAPTIVE_STOPPING_RULE = useAdaptive;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, boolean useGUI, boolean useAdaptive, double targetThreshold) {
        gui = useGUI;
        TARGET_THRESHOLD = targetThreshold;
        USE_TARGET_THRESHOLD = true;
        USE_ADAPTIVE_STOPPING_RULE = useAdaptive;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels, boolean useGUI) {
        gui = useGUI;
        this.varLabels = varLabels;
        haveLabels = true;
        execute(function, draws, burnPeriod, guess);
    }
    
    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels, boolean useGUI, String logFilename) {
        gui = useGUI;
        this.varLabels = varLabels;
        this.logFilename = logFilename;
        haveLabels = true;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels, boolean useGUI, boolean[] ignoreListPassed) {
        gui = useGUI;
        this.varLabels = varLabels;
        haveLabels = true;
        this.ignoreList = ignoreListPassed;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels, boolean useGUI, boolean useAdaptive) {
        gui = useGUI;
        USE_ADAPTIVE_STOPPING_RULE = useAdaptive;
        this.varLabels = varLabels;
        haveLabels = true;
        execute(function, draws, burnPeriod, guess);
    }

    public gibbsLTEGeneralized(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels, boolean useGUI, boolean useAdaptive, double targetThreshold) {
        gui = useGUI;
        USE_ADAPTIVE_STOPPING_RULE = useAdaptive;
        this.varLabels = varLabels;
        haveLabels = true;
        TARGET_THRESHOLD = targetThreshold;
        USE_TARGET_THRESHOLD = true;
        execute(function, draws, burnPeriod, guess);
    }

    private void execute(final mcmcFunction function, int draws, int burnPeriodPassed, double[] guess) {
        if (TUNE_VARIANCE && ACCEPT_ONLY_IMPROVEMENTS) {
            System.out.println("Mutually incompatible parameters: Tuning Variance and Downhill Only");
            System.exit(0);
        }

        if (ignoreList == null) {
            ignoreList = new boolean[guess.length - 1];
            for (int i = 0; i < guess.length - 1; i++) {
                ignoreList[i] = false;
            }
        }

//        for (int i = 0; i < guess.length - 1; i++) {
//            System.out.println("Ignore x[" + (i + 1) + "]: " + ignoreList[i]);
//        }
        long t1 = System.currentTimeMillis();
        updateNumInstances();

        final int burnPeriod = burnPeriodPassed;
        lowestPoint = new double[guess.length];
        updateLowestPoint(guess);
        lowestValue = function.objectiveFunction(guess);

        final int numParams = guess.length;
        point = new double[numParams];
        if (!loadedVar) {
            variance = new double[numParams];
        }
        double[] count = new double[numParams];

        for (int i = 0; i < numParams; i++) {
            point[i] = guess[i];
            if (!loadedVar) {
                variance[i] = Math.max(minimumVariance, Math.abs(point[i] / 30.0));
            }
            count[i] = 0;
        }

        double fold = function.objectiveFunction(point);

        Runnable doSetProgressBarValue = new Runnable() {

            public void run() {
                showParameters = new JCheckBox("Parameter", true);
                showValue = new JCheckBox("Value", true);
                keepGoing = new JCheckBox("Burn-In", false);
                logYAxis = new JCheckBox("Log Y-Axis", false);
                executeNowButton = new JButton("Execute Immediately");
                exitNow = new JButton("Exit");
                shrinkVarianceButton = new JButton("Shrink Variance");
                movingAverageLabel = new JLabel("MA: " + 0);

                chart = ChartFactory.createXYLineChart("Optimization Progress", "Draw", "Objective Function", xyc, PlotOrientation.VERTICAL, true, true, true);

                f = new javax.swing.JFrame("Gibbs LTE Monitor");
                labels = new javax.swing.JLabel[numParams];
                java.text.NumberFormat nf = java.text.NumberFormat.getInstance();
                javax.swing.JPanel panel = new javax.swing.JPanel(true);
                panel.setLayout(new javax.swing.BoxLayout(panel, javax.swing.BoxLayout.Y_AXIS));
                javax.swing.JScrollPane jscr = new javax.swing.JScrollPane(panel);
                javax.swing.JPanel bigPanel = new javax.swing.JPanel();
                bigPanel.setLayout(new java.awt.BorderLayout());
                // bigPanel.add(jscr, java.awt.BorderLayout.CENTER);

                String[] columnNames = {"Parameter", "Value", "MCMC se"};
                String[][] tableData = new String[numParams - 1][3];
                statusTable = new JTable(tableData, columnNames);

                TableModel model = statusTable.getModel();
                bigPanel.add(new JScrollPane(statusTable), BorderLayout.WEST);

                bigPanel.add(new ChartPanel(chart), BorderLayout.EAST);
                javax.swing.JPanel checkPanel = new javax.swing.JPanel();
                checkPanel.setLayout(new java.awt.FlowLayout());
                checkPanel.add(showParameters);
                checkPanel.add(showValue);
                checkPanel.add(keepGoing);
                checkPanel.add(logYAxis);
                checkPanel.add(movingAverageLabel);
                checkPanel.add(executeNowButton);
                checkPanel.add(exitNow);
                checkPanel.add(shrinkVarianceButton);

                bigPanel.add(checkPanel, java.awt.BorderLayout.SOUTH);
                f.getContentPane().add(bigPanel, java.awt.BorderLayout.CENTER);
                bar = new javax.swing.JProgressBar(0, burnPeriod);
                f.getContentPane().add(bar, java.awt.BorderLayout.SOUTH);
                fval = new javax.swing.JLabel("f: unevaluated");
                JPanel topPanel = new JPanel();
                topPanel.add(fval);
                topPanel.add(movingAverageLabel);
                topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
                f.getContentPane().add(topPanel, java.awt.BorderLayout.NORTH);
                for (int j = 1; j < numParams; j++) {
                    labels[j] = new javax.swing.JLabel(nf.format(j - 1)); // names[j-1]+": "+nf.format(point[j])+" (1)");
                    if (haveLabels) {
//                        System.out.println("Labels: "+varLabels.length);
//                        System.out.println("Guess: "+point.length);
//                        System.out.println("Variance: "+variance.length);
                        labels[j].setText(varLabels[j - 1] + ": " + nf.format(point[j]) + " (" + nf.format(variance[j]) + ")");
                    } else {
                        labels[j].setText(nf.format(j - 1) + ": " + nf.format(point[j]) + " (" + nf.format(variance[j]) + ")");
                    }
                    panel.add(labels[j]);
                }
                f.setBounds(100, 50, 500, 400);
                f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                f.pack();
                f.setVisible(true);

                xyc.addSeries(xyBest);
                xyc.addSeries(xyCurrent);

                NumberAxis xAxis = new NumberAxis();
                // xAxis.setAutoRangeStickyZero(false);
                xAxis.setAutoRangeIncludesZero(false);
                chart.getXYPlot().setRangeAxis(0, xAxis);

                // chartDialog.getContentPane().setLayout(new BorderLayout());
                // chartDialog.getContentPane().add(new ChartPanel(chart), BorderLayout.CENTER);
                // chartDialog.setBounds(800, 250, 500, 400);
                // chartDialog.setVisible(true);
                executeNowButton.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent e) {
                        continueExecution = false;
                        executeNowButton.setText("Preparing to Execute...");
                    }
                });
                shrinkVarianceButton.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent e) {
                        for (int i = 0; i < variance.length; i++) {
                            variance[i] = variance[i] / 3.0;
                        }
                    }
                });
                exitNow.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent e) {
                        System.exit(0);
                    }
                });
                logYAxis.addActionListener(new ActionListener() {

                    public void actionPerformed(ActionEvent e) {
                        XYPlot plot = chart.getXYPlot();
                        if (logYAxis.isSelected()) {
                            LogAxis lax = new LogAxis("Ln Objective");
                            lax.setBase(Math.E);
                            plot.setRangeAxis(lax);
                        } else {
                            plot.setRangeAxis(new NumberAxis("Objective"));
                        }
                    }
                });
            }
        };
        if (gui) {
            try {
                SwingUtilities.invokeAndWait(doSetProgressBarValue);
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(0);
            }
        }

        // long t1 = System.currentTimeMillis();
        int ijk = 0;
        fold = function.objectiveFunction(point);
        if (!quiet) {
            System.out.print("Downhill: " + ACCEPT_ONLY_IMPROVEMENTS + " Starting with: " + (-fold) + " ==> ");
            pmUtility.prettyPrint(new Jama.Matrix(guess, 1));
        }
        // System.out.println("burnperiod: "+burnPeriod);
        ProbabilityDistribution normal = new NormalDistribution(0, 1);

        if (gui) {
            updateBarMaximum(burnPeriod);
            for (int j = 1; j < point.length; j++) {
                updateLabel(j, guess[j], variance[j]);
            }
            setBarTextDrawn(true);
        }

        int t = 0;
        int drawsSinceImprovement = 0;
        int m = 0;

        boolean goBurn = true;
        while (goBurn) {
            if (goBurn) {
                // while ((!USE_ADAPTIVE_STOPPING_RULE && ijk < burnPeriod) || (USE_ADAPTIVE_STOPPING_RULE && drawsSinceImprovement <= STOPPING_CRITERION)) {
                if (USE_ADAPTIVE_STOPPING_RULE) {
                    if (drawsSinceImprovement > STOPPING_CRITERION) {
                        goBurn = false;
                    }
                }
                if (USE_TARGET_THRESHOLD) {
                    if (-lowestValue < TARGET_THRESHOLD) {
                        goBurn = false;
                        // System.out.println(lowestValue+" "+TARGET_THRESHOLD);
                    }
                }
                if (ijk > burnPeriod) {
                    goBurn = false;
                }

                double[] pointOld = new double[numParams];
                t++;
                drawsSinceImprovement++;
                // pmUtility.prettyPrint(new Jama.Matrix(point, 1));
                for (int j = 1; j < numParams; j++) {
                    if (!ignoreList[j - 1]) {
                        double pt = point[j] + variance[j] * normal.inverse(Math.random());
                        if (DRAW_ALL_SIMULTANEOUSLY) {
                            for (int lm = 1; lm < numParams; lm++) {
                                pointOld[lm] = point[lm];
                                point[lm] = point[lm] + variance[lm] * normal.inverse(Math.random());
                            }
                        } else {
                            pointOld[j] = point[j];
                            point[j] = pt;
                        }

                        double fnew = function.objectiveFunction(point);
                        // System.out.println("j: " + j + " old: " + (-fold) + " new: " + (-fnew) + " - " + pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1)) + " ");
                        double diff = fnew - fold;
                        double draw = Math.log(Math.random());
                        if (ACCEPT_ONLY_IMPROVEMENTS) {
                            draw = 0;
                        }
                        if (fnew >= lowestValue) {
                            lowestValue = fnew;
                            updateLowestPoint(point);
                            // System.out.println(t + ",x[" + j + "]=" + pt + "," + fnew);
                            drawsSinceImprovement = 0;
                            logUpdate(point, lowestValue);
                        }
                        if (draw <= diff && function.pi(point) > 0.0) {
                            // System.out.println("accept");
                            if (gui) {
                                updateFval(fnew);
                                updateLabel(j, point[j], variance[j]);
                                if (DRAW_ALL_SIMULTANEOUSLY) {
                                    for (int lm = 1; lm < numParams; lm++) {
                                        updateLabel(lm, point[lm], variance[lm]);
                                    }
                                }
                            }
                            fold = fnew;
                        } else {
                            // System.out.println("reject");
                            point[j] = pointOld[j];
                            if (DRAW_ALL_SIMULTANEOUSLY) {
                                for (int lm = 1; lm < numParams; lm++) {
                                    point[lm] = pointOld[lm];
                                }
                            }
                        }
                        if (gui) {
                            if (!Double.isInfinite(transform(-lowestValue))) {
                                xyBest.add(m, transform(-lowestValue));
                            }
                            if (!Double.isInfinite(transform(-fnew))) {
                                xyCurrent.add(m, transform(-fnew));
                            }
                        }
                        m++;
                    }
                }
                ijk++;
                if (gui) {
                    updateBar(ijk);
                    updateBarText((ijk) + " / " + burnPeriod);
                }

            }
        }
        // select variance using the replace one at a time method
        boolean go = true;
        if (!loadedVar) {
            for (int i = 0; i < numParams; i++) {
                variance[i] = Math.max(minimumVariance, Math.abs(point[i] / 30.0));
            }
        } else {
            for (int i = 0; i < numParams; i++) {
                variance[i] = savedVar[i];
            }
        }
        // pmUtility.prettyPrint(new Jama.Matrix(variance, 1));
        /**
         * Tuning the variance of the proposal density.
         */
        while (go && continueExecution && TUNE_VARIANCE) {
            for (int j = 1; j < numParams; j++) {
                // for(int j=3;j<=3;j++) {
                double acceptInside = 0;
                double rejectInside = 0;
                for (int i = 0; i < 25; i++) {
                    double pt = point[j] + variance[j] * normal.inverse(Math.random());
                    double po = point[j];
                    point[j] = pt;
                    double fnew = function.objectiveFunction(point);
                    double diff = fnew - fold;
                    double draw = Math.log(Math.random());
                    if (ACCEPT_ONLY_IMPROVEMENTS) {
                        draw = 0;
                    }
                    if (fnew >= lowestValue) {
                        lowestValue = fnew;
                        updateLowestPoint(point);
                        logUpdate(point, lowestValue);
                    }
                    if (draw <= diff && function.pi(point) > 0.0) {
                        // System.out.println("A--fnew: "+fnew+" diff: "+diff+" draw: "+draw);
                        acceptInside++;
                        if (gui) {
                            if (showValue.isSelected()) {
                                updateFval(fold);
                            }
                            if (showParameters.isSelected()) {
                                updateLabel(j, point[j], variance[j]);
                            }
                        }
                        fold = fnew;
                    } else {
                        // System.out.println("R--fnew: "+fnew+" diff: "+diff+" draw: "+draw);
                        rejectInside++;
                        point[j] = po;
                    }
                    if (gui) {
                        if (!Double.isInfinite(transform(-lowestValue))) {
                            xyBest.add(m, transform(-lowestValue));
                        }
                        if (!Double.isInfinite(transform(-fnew))) {
                            xyCurrent.add(m, transform(-fnew));
                        }
                    }
                    m++;
                }
                double acceptRate = acceptInside / (acceptInside + rejectInside);
                // System.out.print("variance: "+variance[j]+" rate: "+acceptRate+" reject: "+rejectInside+" accept: "+acceptInside+"\n");
                // System.out.print("vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1)));
                double buffer = 0.20;
                // double change = 0.15;
                if (acceptRate > (1.0 + buffer) * acceptance_ratio) {
                    // if(acceptRate>0.45) {
                    variance[j] *= 1.15; // (1.0+change);
                    count[j] = Math.max(0, count[j] - 1);
                }
                if (acceptRate < (1.0 - buffer) * acceptance_ratio) {
                    // if(acceptRate<0.10) {
                    variance[j] *= 0.85; // (1.0-change);
                    count[j] = Math.max(0, count[j] - 1);
                }
                if (acceptRate == 0) {
                    variance[j] *= 0.5;
                }
                if (acceptRate == 1) {
                    variance[j] *= 2;
                }
                if (acceptRate <= acceptance_ratio * (1.0 + buffer) && acceptRate >= acceptance_ratio * (1.0 - buffer)) {
                    // if(acceptRate<0.45 && acceptRate>0.10) {
                    count[j]++;
                }
                if (gui) {
                    if (showParameters.isSelected()) {
                        updateLabel(j, point[j], variance[j], acceptRate);
                    }
                }
            }
            go = false;
            for (int i = 1; i < numParams; i++) {
                if (!go) {
                    if (count[i] >= 3) {
                        go = false;
                    } else {
                        go = true;
                    }
                }
            }
        }
        // System.out.print("Post variance selection: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1))+"\n");
        // pmUtility.prettyPrint(new Jama.Matrix(variance, 1));
        if (gui) {
            updateBar(0);
            changeNowButtonText("Executing...");
        }
        // now run chain
        chain = new Jama.Matrix(draws, numParams);
        int accept = 0;
        int reject = 0;
        if (gui) {
            updateBarMaximum(draws);
            setBarTextDrawn(true);
        }
        // for (int i = 0; i < draws; i++) {
        int drawsSinceLastImprovement = 0;
        int i = -1;
        boolean drawingOfTheThree = true;
        while (drawingOfTheThree) {
            if (drawingOfTheThree) {
                i++;
                if (USE_TARGET_THRESHOLD) {
                    if (-lowestValue < TARGET_THRESHOLD) {
                        drawingOfTheThree = false;
                        // System.out.println(lowestValue+" "+TARGET_THRESHOLD);
                    }
                } else {
                    if (USE_ADAPTIVE_STOPPING_RULE && drawsSinceLastImprovement > STOPPING_CRITERION) {
                        drawingOfTheThree = false;
                    }
                    if (!USE_ADAPTIVE_STOPPING_RULE && i >= draws - 1) {
                        drawingOfTheThree = false;
                    }
                }

                t++;
                if (gui) {
                    if (!USE_ADAPTIVE_STOPPING_RULE) {
                        updateBar(i + 1);
                        updateBarText((i + 1) + " / " + draws);
                    }
                }
                drawsSinceLastImprovement++;
                double[] pointOld = new double[point.length];
                for (int j = 1; j < numParams; j++) {
                    if (!ignoreList[j - 1]) {
                        // double pt = point[j]+cauchy.inverse(Math.random());

                        double pt = point[j] + variance[j] * normal.inverse(Math.random());
                        double po = point[j];

                        if (DRAW_ALL_SIMULTANEOUSLY) {
                            for (int lm = 1; lm < numParams; lm++) {
                                pointOld[lm] = point[lm];
                                point[lm] = point[lm] + variance[lm] * normal.inverse(Math.random());
                            }
                        }

                        point[j] = pt;
                        double fnew = function.objectiveFunction(point);
                        double diff = fnew - fold;
                        double draw = Math.log(Math.random());
                        if (ACCEPT_ONLY_IMPROVEMENTS) {
                            draw = 0;
                        }

//                    if (GENERATE_SET_SAMPLES) {
//                        for (int k = 0; k < numParams; k++) {
//                            chain.set(i, k, point[k]);
//                        }
//                    }
                        if (draw <= diff && function.pi(point) > 0.0) {
                            accept++;
                            fold = fnew;
                            if (gui) {
                                if (showValue.isSelected()) {
                                    updateFval(fold);
                                }
                                if (showParameters.isSelected()) {
                                    updateLabel(j, point[j], variance[j]);
                                    if (DRAW_ALL_SIMULTANEOUSLY) {
                                        for (int lm = 1; lm < numParams; lm++) {
                                            updateLabel(lm, point[lm], variance[lm]);
                                        }
                                    }
                                }
                            }
                            // System.out.print("f: "+fnew+" vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1)));
                        } else {
                            reject++;
                            point[j] = po;
                            if (DRAW_ALL_SIMULTANEOUSLY) {
                                for (int lm = 1; lm < numParams; lm++) {
                                    point[lm] = pointOld[lm];
                                }
                            }
                        }

                        if (fnew > lowestValue) {
                            lowestValue = fnew;
                            updateLowestPoint(point);
                            drawsSinceLastImprovement = 0;
                            logUpdate(point, lowestValue);
                            // System.out.println(t + ",x[" + j + "]=" + pt + "," + fnew);
                        }

                        if (gui) {
                            if (!Double.isInfinite(transform(-lowestValue))) {
                                xyBest.add(m, transform(-lowestValue));
                            }
                            if (!Double.isInfinite(transform(-fnew))) {
                                xyCurrent.add(m, transform(-fnew));
                            }
                        }
                        m++;
                    }
                }
                if (!USE_ADAPTIVE_STOPPING_RULE && !USE_TARGET_THRESHOLD) {
                    for (int k = 0; k < numParams; k++) {
                        chain.set(i, k, point[k]);
                    }
                }
            }
        }
        // System.out.print("Final vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1))+"\n");
        // System.out.print("Reject: "+reject+" accept: "+accept+"\n");
        if (gui) {
            Runnable updateLabels = new Runnable() {

                public void run() {
                    f.dispose();
                }
            };
            SwingUtilities.invokeLater(updateLabels);
        }
        long t2 = System.currentTimeMillis();
        if (!quiet) {
            System.out.print((t2 - t1) + " ms t: " + t + " Lowest point stored: " + (-lowestValue) + " ");
            pmUtility.prettyPrint(new Jama.Matrix(lowestPoint, 1));
        }
    }

    public void updateLowestPoint(double[] x) {
        for (int i = 0; i < x.length; i++) {
            lowestPoint[i] = x[i];
        }
    }

    public void updateFval(final double x) {
        if (gui) {
            Runnable updateLabels = new Runnable() {

                public void run() {
                    fval.setText(String.format("f: %,g", (-x)));
                }
            };
            SwingUtilities.invokeLater(updateLabels);
        }
    }

    private void logUpdate(double[] bestX, double fValueBest) {
        if (logFilename != null) {
            try {
                BufferedWriter writer = new BufferedWriter(new FileWriter(logFilename, true));
                writer.append(-fValueBest + ": ");
                writer.append(pmUtility.stringPrettyPrint(new Jama.Matrix(bestX, 1)) + "\n");
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(0);
            }
        }
    }

    private void changeNowButtonText(final String s) {
        if (gui) {
            Runnable updateLabels = new Runnable() {

                public void run() {
                    executeNowButton.setText(s);
                }
            };
            SwingUtilities.invokeLater(updateLabels);
        }
    }

    public void updateLabel(final int j, final double point, final double variance) {
        if (gui) {
            Runnable updateLabels = new Runnable() {

                public void run() {
                    if (haveLabels) {
                        labels[j].setText(varLabels[j - 1] + ": " + nf.format(point) + " (" + nf.format(variance) + ")");
                        statusTable.getModel().setValueAt(varLabels[j - 1], j - 1, 0);
                    } else {
                        labels[j].setText(nf.format(j - 1) + ": " + nf.format(point) + " (" + nf.format(variance) + ")");
                        statusTable.getModel().setValueAt("x[" + j + "]:", j - 1, 0);
                    }
                    statusTable.getModel().setValueAt(nf.format(point), j - 1, 1);
                    statusTable.getModel().setValueAt(nf.format(variance), j - 1, 2);
                }
            };
            SwingUtilities.invokeLater(updateLabels);
        }
    }

    private void updateBar(final int i) {
        if (gui) {
            Runnable updateBar = new Runnable() {

                public void run() {
                    bar.setValue(i);
                }
            };
            SwingUtilities.invokeLater(updateBar);
        }
    }

    public void notifyBurnProgessed() {
        burnProgress++;
        // System.out.println(burnProgress);
        updateBar(burnProgress);
    }

    private void updateBarMaximum(final int i) {
        if (gui) {
            Runnable updateBar = new Runnable() {

                public void run() {
                    bar.setMaximum(i);
                }
            };
            SwingUtilities.invokeLater(updateBar);
        }
    }

    public void updateLabel(final int j, final double point, final double variance, final double acceptRate) {
        if (gui) {
            Runnable updateLabels = new Runnable() {

                public void run() {
                    if (haveLabels) {
                        labels[j].setText(varLabels[j - 1] + ": " + nf.format(point) + " (" + nf.format(variance) + ") rate: " + nf.format(acceptRate));
                        statusTable.getModel().setValueAt(varLabels[j - 1], j - 1, 0);
                    } else {
                        labels[j].setText(nf.format(j - 1) + ": " + nf.format(point) + " (" + nf.format(variance) + ") rate: " + nf.format(acceptRate));
                        statusTable.getModel().setValueAt("x[" + j + "]:", j - 1, 0);
                    }
                    statusTable.getModel().setValueAt(nf.format(point) + " rate: " + nf.format(acceptRate), j - 1, 1);
                    statusTable.getModel().setValueAt(nf.format(variance), j - 1, 2);
                }
            };
            SwingUtilities.invokeLater(updateLabels);
        }
    }

    public Jama.Matrix getChain() {
        return chain;
    }

    public double[] getLowestPoint() {
        return lowestPoint;
    }

    private void updateBarText(final String string) {
        if (gui) {
            Runnable task = new Runnable() {

                public void run() {
                    bar.setString(string);
                }
            };
            SwingUtilities.invokeLater(task);
        }
    }

    private void setBarTextDrawn(final boolean b) {
        if (gui) {
            Runnable task = new Runnable() {

                public void run() {
                    bar.setStringPainted(b);
                }
            };
            SwingUtilities.invokeLater(task);
        }
    }

    private synchronized void updateNumInstances() {
        Runnable task = new Runnable() {

            public void run() {
                NUM_INSTANCES++;
                if (!quiet) {
                    System.out.println("Ding! (" + NUM_INSTANCES + ")");
                }
            }
        };
        SwingUtilities.invokeLater(task);
    }

    private double transform(double x) {
        boolean logTransformation = false;
        if (logTransformation) {
            if (x > 0) {
                return Math.log(1 + x);
            } else {
                return x;
            }
        } else {
            return x;
        }
    }
}
