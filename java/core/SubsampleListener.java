package core;

import java.util.ArrayList;

/**
 * Callback interface for receiving real-time subsampling progress from TreeMoment.
 * Allows the GUI to display the distribution building up as subsamples complete.
 */
public interface SubsampleListener {

    /**
     * Called when subsampling begins for a parameter. Returns a token identifying
     * this subsampling run so incremental updates can target the right chart.
     *
     * @param treeIndex the index of the tree in the forest
     * @param parameterIndex the parameter being tested (k)
     * @param Tn the observed test statistic on the full sample
     * @param totalSubsamples total number of subsamples planned
     * @return a token object to pass to onSubsampleProgress and onSubsampleComplete
     */
    Object onSubsampleStart(int treeIndex, int parameterIndex, double Tn, int totalSubsamples);

    /**
     * Called periodically as subsamples complete, to update the live chart.
     *
     * @param token the token returned from onSubsampleStart
     * @param stats all subsampled test statistics so far
     * @param completedCount how many subsamples have completed
     */
    void onSubsampleProgress(Object token, ArrayList<Double> stats, int completedCount);

    /**
     * Called when subsampling for one parameter is fully complete.
     *
     * @param token the token returned from onSubsampleStart
     * @param stats the final subsampled test statistics
     * @param Tn the observed test statistic on the full sample
     * @param criticalValue the 95th percentile critical value
     * @param pValue the p-value
     */
    void onSubsampleComplete(Object token, ArrayList<Double> stats, double Tn, double criticalValue, double pValue);
}
