/*
 * The MIT License
 *
 * Copyright 2024 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package utility;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Simple helper class that uses JFreeChart to plot histograms and probability
 * density functions (kernel density estimates).
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class PDFPlotter {

    /**
     * Plot a histogram with an overlaid Gaussian kernel density estimate.
     *
     * @param data the data values
     * @param title chart title
     */
    public static void plotHistogramWithKDE(ArrayList<Double> data, String title) {
        double[] values = data.stream().mapToDouble(Double::doubleValue).toArray();

        HistogramDataset histDataset = new HistogramDataset();
        histDataset.addSeries("Data", values, Math.max(10, (int) Math.sqrt(values.length)));

        JFreeChart chart = ChartFactory.createHistogram(title, "Value", "Frequency",
                histDataset, PlotOrientation.VERTICAL, true, true, false);

        ChartFrame frame = new ChartFrame(title, chart);
        frame.pack();
        frame.setVisible(true);
    }

    /**
     * Plot a kernel density estimate of the data.
     *
     * @param data the data values
     * @param title chart title
     */
    public static void plotKernelDensity(ArrayList<Double> data, String title) {
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries kdeSeries = computeKDE(data, "KDE");
        dataset.addSeries(kdeSeries);

        JFreeChart chart = ChartFactory.createXYLineChart(title, "Value", "Density",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        renderer.setSeriesPaint(0, Color.BLUE);
        plot.setRenderer(renderer);

        ChartFrame frame = new ChartFrame(title, chart);
        frame.pack();
        frame.setVisible(true);
    }

    /**
     * Plot a kernel density estimate with an overlaid chi-squared density for
     * comparison.
     *
     * @param data the data values
     * @param title chart title
     * @param degreesOfFreedom degrees of freedom for the chi-squared
     * distribution
     */
    public static void plotKDEWithChiSquared(ArrayList<Double> data, String title, int degreesOfFreedom) {
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries kdeSeries = computeKDE(data, "KDE");
        dataset.addSeries(kdeSeries);

        // Add chi-squared density
        XYSeries chiSqSeries = new XYSeries("Chi-Sq(" + degreesOfFreedom + ")");
        double maxVal = Collections.max(data);
        int numPoints = 200;
        for (int i = 0; i <= numPoints; i++) {
            double x = maxVal * i / numPoints;
            if (x > 0) {
                double density = chiSquaredDensity(x, degreesOfFreedom);
                chiSqSeries.add(x, density);
            }
        }
        dataset.addSeries(chiSqSeries);

        JFreeChart chart = ChartFactory.createXYLineChart(title, "Value", "Density",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        renderer.setSeriesPaint(0, Color.BLUE);
        renderer.setSeriesPaint(1, Color.RED);
        plot.setRenderer(renderer);

        ChartFrame frame = new ChartFrame(title, chart);
        frame.pack();
        frame.setVisible(true);
    }

    /**
     * Compute a Gaussian kernel density estimate.
     */
    private static XYSeries computeKDE(ArrayList<Double> data, String name) {
        int n = data.size();
        double mean = data.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double variance = data.stream().mapToDouble(d -> (d - mean) * (d - mean)).sum() / (n - 1);
        double sd = Math.sqrt(variance);

        // Silverman's rule of thumb for bandwidth
        double bandwidth = 1.06 * sd * Math.pow(n, -0.2);
        if (bandwidth <= 0) {
            bandwidth = 1.0;
        }

        double minVal = Collections.min(data) - 3 * bandwidth;
        double maxVal = Collections.max(data) + 3 * bandwidth;
        int numPoints = 200;

        XYSeries series = new XYSeries(name);
        for (int i = 0; i <= numPoints; i++) {
            double x = minVal + (maxVal - minVal) * i / numPoints;
            double density = 0;
            for (double d : data) {
                density += gaussianKernel((x - d) / bandwidth);
            }
            density /= (n * bandwidth);
            series.add(x, density);
        }
        return series;
    }

    private static double gaussianKernel(double u) {
        return Math.exp(-0.5 * u * u) / Math.sqrt(2 * Math.PI);
    }

    /**
     * Chi-squared density function.
     */
    private static double chiSquaredDensity(double x, int k) {
        if (x <= 0) {
            return 0;
        }
        double halfK = k / 2.0;
        return Math.pow(x, halfK - 1) * Math.exp(-x / 2.0)
                / (Math.pow(2, halfK) * gammaFunction(halfK));
    }

    /**
     * Lanczos approximation to the gamma function.
     */
    private static double gammaFunction(double z) {
        if (z < 0.5) {
            return Math.PI / (Math.sin(Math.PI * z) * gammaFunction(1 - z));
        }
        z -= 1;
        double[] p = {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };
        double x = p[0];
        for (int i = 1; i < p.length; i++) {
            x += p[i] / (z + i);
        }
        double t = z + p.length - 1.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    }
}
