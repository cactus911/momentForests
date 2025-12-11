package core;

import JSci.maths.statistics.ChiSqrDistribution;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

public class PDFPlotter extends JFrame {
    
    /**
     * Plot histogram (empirical PDF) from ArrayList of Doubles
     */
    public static void plotHistogram(ArrayList<Double> data, String title) {
        // Convert ArrayList to double array
        double[] dataArray = data.stream()
            .mapToDouble(Double::doubleValue)
            .toArray();
        
        // Calculate optimal number of bins using Freedman-Diaconis rule
        int bins = calculateOptimalBins(dataArray);
        
        // Create histogram dataset
        HistogramDataset dataset = new HistogramDataset();
        dataset.addSeries("PDF", dataArray, bins);
        
        // Create chart
        JFreeChart chart = ChartFactory.createHistogram(
            title,
            "Value",
            "Probability Density",
            dataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        // Customize chart appearance
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
    
    /**
     * Plot smooth PDF using kernel density estimation
     */
    public static void plotKernelDensity(ArrayList<Double> data, String title) {
        double[] dataArray = data.stream()
            .mapToDouble(Double::doubleValue)
            .toArray();
        
        // Calculate KDE
        XYSeries series = new XYSeries("PDF");
        double[] kde = kernelDensityEstimation(dataArray, 200);
        double min = Arrays.stream(dataArray).min().orElse(0);
        double max = Arrays.stream(dataArray).max().orElse(0);
        double step = (max - min) / (kde.length - 1);
        
        for (int i = 0; i < kde.length; i++) {
            double x = min + i * step;
            series.add(x, kde[i]);
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        
        // Create chart
        JFreeChart chart = ChartFactory.createXYLineChart(
            title,
            "Value",
            "Probability Density",
            dataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        // Customize appearance
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesPaint(0, new Color(0, 100, 200));
        renderer.setSeriesStroke(0, new BasicStroke(2.0f));
        plot.setRenderer(renderer);
        
        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
    
    /**
     * Plot both histogram and KDE overlaid
     */
    public static void plotHistogramWithKDE(ArrayList<Double> data, String title) {
        double[] dataArray = data.stream()
            .mapToDouble(Double::doubleValue)
            .toArray();
        
        int bins = calculateOptimalBins(dataArray);
        
        // Create histogram dataset
        HistogramDataset histDataset = new HistogramDataset();
        histDataset.addSeries("Histogram", dataArray, bins);
        
        // Create KDE dataset
        XYSeries kdeSeries = new XYSeries("KDE");
        double[] kde = kernelDensityEstimation(dataArray, 200);
        double min = Arrays.stream(dataArray).min().orElse(0);
        double max = Arrays.stream(dataArray).max().orElse(0);
        double step = (max - min) / (kde.length - 1);
        
        // Scale KDE to match histogram area
        double binWidth = (max - min) / bins;
        double scaleFactor = dataArray.length * binWidth;
        
        for (int i = 0; i < kde.length; i++) {
            double x = min + i * step;
            kdeSeries.add(x, kde[i] * scaleFactor);
        }
        
        XYSeriesCollection kdeDataset = new XYSeriesCollection(kdeSeries);
        
        // Create chart with histogram
        JFreeChart chart = ChartFactory.createHistogram(
            title,
            "Value",
            "Frequency",
            histDataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        XYPlot plot = chart.getXYPlot();
        plot.setDataset(1, kdeDataset);
        
        // Customize renderers
        XYLineAndShapeRenderer kdeRenderer = new XYLineAndShapeRenderer();
        kdeRenderer.setSeriesLinesVisible(0, true);
        kdeRenderer.setSeriesShapesVisible(0, false);
        kdeRenderer.setSeriesPaint(0, Color.RED);
        kdeRenderer.setSeriesStroke(0, new BasicStroke(2.5f));
        plot.setRenderer(1, kdeRenderer);
        
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
    
    // Helper: Kernel Density Estimation using Gaussian kernel
    private static double[] kernelDensityEstimation(double[] data, int points) {
        double min = Arrays.stream(data).min().orElse(0);
        double max = Arrays.stream(data).max().orElse(0);
        double bandwidth = calculateBandwidth(data);
        
        double[] kde = new double[points];
        double step = (max - min) / (points - 1);
        
        for (int i = 0; i < points; i++) {
            double x = min + i * step;
            double sum = 0;
            
            for (double datum : data) {
                double u = (x - datum) / bandwidth;
                sum += gaussianKernel(u);
            }
            
            kde[i] = sum / (data.length * bandwidth);
        }
        
        return kde;
    }
    
    private static double gaussianKernel(double u) {
        return Math.exp(-0.5 * u * u) / Math.sqrt(2 * Math.PI);
    }
    
    private static double calculateBandwidth(double[] data) {
        double stddev = standardDeviation(data);
        double n = data.length;
        // Silverman's rule of thumb
        return 1.06 * stddev * Math.pow(n, -0.2);
    }
    
    private static double standardDeviation(double[] data) {
        double mean = Arrays.stream(data).average().orElse(0);
        double variance = Arrays.stream(data)
            .map(x -> Math.pow(x - mean, 2))
            .average()
            .orElse(0);
        return Math.sqrt(variance);
    }
    
    private static int calculateOptimalBins(double[] data) {
        // Freedman-Diaconis rule
        double[] sorted = Arrays.copyOf(data, data.length);
        Arrays.sort(sorted);
        
        double q1 = percentile(sorted, 25);
        double q3 = percentile(sorted, 75);
        double iqr = q3 - q1;
        
        double n = data.length;
        double h = 2 * iqr / Math.pow(n, 1.0/3.0);
        double range = sorted[sorted.length - 1] - sorted[0];
        
        return Math.max(10, (int) Math.ceil(range / h));
    }
    
    private static double percentile(double[] sortedData, double percentile) {
        int n = sortedData.length;
        double index = percentile / 100.0 * (n - 1);
        int lower = (int) Math.floor(index);
        int upper = (int) Math.ceil(index);
        
        if (lower == upper) {
            return sortedData[lower];
        }
        
        double weight = index - lower;
        return sortedData[lower] * (1 - weight) + sortedData[upper] * weight;
    }
    
    /**
     * Plot histogram with chi-squared distribution overlay
     */
    public static void plotHistogramWithChiSquared(ArrayList<Double> data, String title, int degreesOfFreedom) {
        double[] dataArray = data.stream()
            .mapToDouble(Double::doubleValue)
            .toArray();
        
        int bins = calculateOptimalBins(dataArray);
        
        // Create histogram dataset
        HistogramDataset histDataset = new HistogramDataset();
        histDataset.addSeries("Empirical Data", dataArray, bins);
        
        // Create chi-squared theoretical distribution
        XYSeries chiSqSeries = new XYSeries("χ²(" + degreesOfFreedom + ")");
        double min = Math.max(0, Arrays.stream(dataArray).min().orElse(0));
        double max = Arrays.stream(dataArray).max().orElse(0);
        int points = 500;
        double step = (max - min) / (points - 1);
        
        // Scale chi-squared to match histogram area
        double binWidth = (max - min) / bins;
        double scaleFactor = dataArray.length * binWidth;
        
        JSci.maths.statistics.ChiSqrDistribution chi = new ChiSqrDistribution(degreesOfFreedom);
        
        for (int i = 0; i < points; i++) {
            double x = min + i * step;
            double chiSqValue = chi.probability(x);
            chiSqSeries.add(x, chiSqValue * scaleFactor);
        }
        
        XYSeriesCollection chiSqDataset = new XYSeriesCollection(chiSqSeries);
        
        // Create chart
        JFreeChart chart = ChartFactory.createHistogram(
            title,
            "Value",
            "Frequency",
            histDataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        XYPlot plot = chart.getXYPlot();
        plot.setDataset(1, chiSqDataset);
        
        // Customize renderers
        XYLineAndShapeRenderer chiSqRenderer = new XYLineAndShapeRenderer();
        chiSqRenderer.setSeriesLinesVisible(0, true);
        chiSqRenderer.setSeriesShapesVisible(0, false);
        chiSqRenderer.setSeriesPaint(0, new Color(200, 0, 0));
        chiSqRenderer.setSeriesStroke(0, new BasicStroke(2.5f));
        plot.setRenderer(1, chiSqRenderer);
        
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
    
    /**
     * Plot KDE with chi-squared distribution overlay
     */
    public static void plotKDEWithChiSquared(ArrayList<Double> data, String title, int degreesOfFreedom) {
        double[] dataArray = data.stream()
            .mapToDouble(Double::doubleValue)
            .toArray();
        
        // Calculate KDE
        XYSeries kdeSeries = new XYSeries("Empirical KDE");
        double[] kde = kernelDensityEstimation(dataArray, 200);
        double min = Math.max(0, Arrays.stream(dataArray).min().orElse(0));
        double max = Arrays.stream(dataArray).max().orElse(0);
        double step = (max - min) / (kde.length - 1);
        
        for (int i = 0; i < kde.length; i++) {
            double x = min + i * step;
            kdeSeries.add(x, kde[i]);
        }
        
        // Create chi-squared series
        XYSeries chiSqSeries = new XYSeries("χ²(" + degreesOfFreedom + ")");
        int points = 500;
        step = (max - min) / (points - 1);
        
        for (int i = 0; i < points; i++) {
            double x = min + i * step;
            double chiSqValue = chiSquaredPDF(x, degreesOfFreedom);
            chiSqSeries.add(x, chiSqValue);
        }
        
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(kdeSeries);
        dataset.addSeries(chiSqSeries);
        
        // Create chart
        JFreeChart chart = ChartFactory.createXYLineChart(
            title,
            "Value",
            "Probability Density",
            dataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        // Customize appearance
        XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(0, true);
        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesPaint(0, new Color(0, 100, 200));
        renderer.setSeriesStroke(0, new BasicStroke(2.0f));
        
        renderer.setSeriesLinesVisible(1, true);
        renderer.setSeriesShapesVisible(1, false);
        renderer.setSeriesPaint(1, new Color(200, 0, 0));
        renderer.setSeriesStroke(1, new BasicStroke(2.5f, BasicStroke.CAP_BUTT, 
            BasicStroke.JOIN_MITER, 10.0f, new float[]{10.0f}, 0.0f));
        
        plot.setRenderer(renderer);
        
        // Display chart
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    }
    
    /**
     * Chi-squared PDF: f(x; k) = (1/(2^(k/2) * Γ(k/2))) * x^(k/2-1) * e^(-x/2)
     */
    private static double chiSquaredPDF(double x, int k) {
        if (x <= 0) return 0;
        
        double kHalf = k / 2.0;
        double coefficient = 1.0 / (Math.pow(2, kHalf) * gamma(kHalf));
        double power = Math.pow(x, kHalf - 1);
        double exponential = Math.exp(-x / 2.0);
        
        return coefficient * power * exponential;
    }
    
    /**
     * Gamma function approximation using Stirling's formula
     */
    private static double gamma(double z) {
        if (z == 1.0) return 1.0;
        if (z == 0.5) return Math.sqrt(Math.PI);
        if (z == 1.5) return Math.sqrt(Math.PI) / 2.0;
        if (z == 2.0) return 1.0;
        if (z == 2.5) return 1.5 * Math.sqrt(Math.PI) / 2.0;
        
        // Lanczos approximation for other values
        double[] coef = {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };
        
        if (z < 0.5) {
            return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
        }
        
        z -= 1;
        double x = coef[0];
        for (int i = 1; i < coef.length; i++) {
            x += coef[i] / (z + i);
        }
        
        double t = z + coef.length - 1.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    }
    
    public static void main(String[] args) {
        // Generate sample data from chi-squared distribution
        Random rand = new Random(42);
        ArrayList<Double> chiSqData = new ArrayList<>();
        int df = 5; // degrees of freedom
        
        // Generate chi-squared samples (sum of squared normals)
        for (int i = 0; i < 1000; i++) {
            double sum = 0;
            for (int j = 0; j < df; j++) {
                double z = rand.nextGaussian();
                sum += z * z;
            }
            chiSqData.add(sum);
        }
        
        // Plot histogram with chi-squared overlay
        SwingUtilities.invokeLater(() -> {
            plotHistogramWithChiSquared(chiSqData, 
                "Chi-Squared Distribution (df=" + df + ") - Histogram", df);
        });
        
        // Plot KDE with chi-squared overlay
        SwingUtilities.invokeLater(() -> {
            plotKDEWithChiSquared(chiSqData, 
                "Chi-Squared Distribution (df=" + df + ") - KDE", df);
        });
        
        // Original examples with normal distribution
        ArrayList<Double> normalData = new ArrayList<>();
        for (int i = 0; i < 1000; i++) {
            normalData.add(rand.nextGaussian() * 15 + 50);
        }
        
        SwingUtilities.invokeLater(() -> {
            plotHistogram(normalData, "Normal Distribution - Histogram");
        });
    }
}