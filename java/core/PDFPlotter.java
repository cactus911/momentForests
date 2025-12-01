package core;

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
    
    public static void main(String[] args) {
        // Generate sample data (normal distribution)
        Random rand = new Random(42);
        ArrayList<Double> data = new ArrayList<>();
        for (int i = 0; i < 1000; i++) {
            data.add(rand.nextGaussian() * 15 + 50);
        }
        
        // Plot histogram
        SwingUtilities.invokeLater(() -> {
            plotHistogram(data, "Probability Distribution - Histogram");
        });
        
        // Plot smooth KDE
        SwingUtilities.invokeLater(() -> {
            plotKernelDensity(data, "Probability Distribution - Kernel Density");
        });
        
        // Plot both overlaid
        SwingUtilities.invokeLater(() -> {
            plotHistogramWithKDE(data, "PDF - Histogram with KDE Overlay");
        });
    }
}