package core;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * A single scrollable panel that shows a flowing record of all subsampling
 * activity: text labels, live-updating KDE density charts, and result labels,
 * interleaved top to bottom. Nothing is ever cleared.
 */
public class SubsampleVisualizationPanel extends JPanel implements SubsampleListener {

    private final JPanel flowPanel;
    private final JScrollPane scrollPane;

    private static final int CHART_HEIGHT = 120;
    private static final int CHART_WIDTH = 350;
    private volatile int currentMonteCarloIndex = -1;

    public SubsampleVisualizationPanel() {
        setLayout(new BorderLayout());

        flowPanel = new JPanel();
        flowPanel.setLayout(new BoxLayout(flowPanel, BoxLayout.Y_AXIS));
        scrollPane = new JScrollPane(flowPanel);
        scrollPane.getVerticalScrollBar().setUnitIncrement(32);

        add(scrollPane, BorderLayout.CENTER);
    }

    /**
     * Returns a dummy JTextArea for compatibility with code that calls
     * vizPanel.getLogArea(). Anything appended here gets added as a label
     * in the flow panel instead.
     */
    public JTextArea getLogArea() {
        // Return a proxy JTextArea that redirects appends to the flow panel
        return new FlowProxyTextArea(this);
    }

    /** Set the current Monte Carlo iteration index for labeling charts. */
    public void setMonteCarloIndex(int mcIndex) {
        this.currentMonteCarloIndex = mcIndex;
    }

    /** Add a text line into the flowing scroll. */
    public void appendText(String text) {
        SwingUtilities.invokeLater(() -> {
            boolean wasBottom = isScrolledToBottom();
            JLabel lbl = new JLabel("  " + text.trim());
            lbl.setFont(new java.awt.Font("Monospaced", java.awt.Font.PLAIN, 11));
            lbl.setAlignmentX(LEFT_ALIGNMENT);
            lbl.setMaximumSize(new Dimension(Integer.MAX_VALUE, 16));
            lbl.setBorder(BorderFactory.createEmptyBorder(1, 4, 1, 4));
            flowPanel.add(lbl);
            flowPanel.revalidate();
            if (wasBottom) smartAutoScroll();
        });
    }

    private boolean isScrolledToBottom() {
        JScrollBar vbar = scrollPane.getVerticalScrollBar();
        int extent = vbar.getModel().getExtent();
        return (vbar.getValue() + extent >= vbar.getMaximum() - 40);
    }

    private void smartAutoScroll() {
        SwingUtilities.invokeLater(() -> {
            JScrollBar vbar = scrollPane.getVerticalScrollBar();
            vbar.setValue(vbar.getMaximum());
        });
    }

    // ---- Token ----

    private static class SubsampleToken {
        final int treeIndex;
        final int parameterIndex;
        final double Tn;
        final int totalSubsamples;
        final int mcIndex;
        ChartPanel chartPanel;
        JFreeChart chart;
        JPanel groupPanel;
        JLabel headerLabel;
        JProgressBar progressBar;

        SubsampleToken(int treeIndex, int parameterIndex, double Tn, int totalSubsamples, int mcIndex) {
            this.treeIndex = treeIndex;
            this.parameterIndex = parameterIndex;
            this.Tn = Tn;
            this.totalSubsamples = totalSubsamples;
            this.mcIndex = mcIndex;
        }

        String prefix() {
            return mcIndex >= 0 ? "MC #" + mcIndex + " | " : "";
        }
    }

    // ---- Listener methods ----

    @Override
    public Object onSubsampleStart(int treeIndex, int parameterIndex, double Tn, int totalSubsamples) {
        SubsampleToken token = new SubsampleToken(treeIndex, parameterIndex, Tn, totalSubsamples, currentMonteCarloIndex);
        boolean wasBottom = isScrolledToBottom();

        SwingUtilities.invokeLater(() -> {
            // Empty XY chart (no title — header label serves that role)
            XYSeriesCollection dataset = new XYSeriesCollection();
            JFreeChart chart = ChartFactory.createXYLineChart(
                    null, null, null,
                    dataset, PlotOrientation.VERTICAL, false, false, false);

            // Anti-alias everything
            chart.setAntiAlias(true);
            chart.setTextAntiAlias(true);

            XYPlot plot = chart.getXYPlot();
            plot.setBackgroundPaint(Color.WHITE);
            plot.setDomainGridlinePaint(new Color(230, 230, 230));
            plot.setRangeGridlinePaint(new Color(230, 230, 230));
            plot.getDomainAxis().setTickLabelFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 9));
            plot.getRangeAxis().setTickLabelFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 9));

            // Tn marker
            ValueMarker tnMarker = new ValueMarker(Tn);
            tnMarker.setPaint(Color.BLACK);
            tnMarker.setStroke(new BasicStroke(1.5f, BasicStroke.CAP_BUTT,
                    BasicStroke.JOIN_MITER, 10.0f, new float[]{5.0f}, 0.0f));
            plot.addDomainMarker(tnMarker);

            ChartPanel cp = new ChartPanel(chart);
            cp.setPreferredSize(new Dimension(CHART_WIDTH, CHART_HEIGHT));
            cp.setMaximumSize(new Dimension(CHART_WIDTH, CHART_HEIGHT));
            cp.setMinimumSize(new Dimension(200, CHART_HEIGHT));

            token.chart = chart;
            token.chartPanel = cp;

            // Group: header + chart (result added on complete)
            JPanel group = new JPanel();
            group.setLayout(new BoxLayout(group, BoxLayout.Y_AXIS));
            group.setBorder(BorderFactory.createCompoundBorder(
                    BorderFactory.createMatteBorder(0, 0, 1, 0, new Color(210, 210, 210)),
                    BorderFactory.createEmptyBorder(2, 2, 2, 2)));

            JLabel header = new JLabel("  " + token.prefix() + "Tree " + treeIndex + ", k=" + parameterIndex
                    + " \u2014 0/" + totalSubsamples);
            header.setFont(new java.awt.Font("SansSerif", java.awt.Font.BOLD, 10));
            header.setAlignmentX(LEFT_ALIGNMENT);
            header.setMaximumSize(new Dimension(Integer.MAX_VALUE, 18));
            token.headerLabel = header;

            cp.setAlignmentX(LEFT_ALIGNMENT);

            JProgressBar pbar = new JProgressBar(0, totalSubsamples);
            pbar.setValue(0);
            pbar.setStringPainted(true);
            pbar.setString("0/" + totalSubsamples);
            pbar.setFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 9));
            pbar.setAlignmentX(LEFT_ALIGNMENT);
            pbar.setMaximumSize(new Dimension(CHART_WIDTH, 14));
            pbar.setPreferredSize(new Dimension(CHART_WIDTH, 14));
            token.progressBar = pbar;

            group.add(header);
            group.add(cp);
            group.add(pbar);
            group.setMaximumSize(new Dimension(Integer.MAX_VALUE, CHART_HEIGHT + 55));

            token.groupPanel = group;
            flowPanel.add(group);
            flowPanel.revalidate();
            flowPanel.repaint();

            if (wasBottom) smartAutoScroll();
        });

        return token;
    }

    @Override
    public void onSubsampleProgress(Object tokenObj, ArrayList<Double> stats, int completedCount) {
        SubsampleToken token = (SubsampleToken) tokenObj;
        if (token.chart == null || stats.size() < 15) return;

        ArrayList<Double> copy = new ArrayList<>(stats);
        SwingUtilities.invokeLater(() -> {
            updateDensityPlot(token, copy, completedCount, false, 0, 0);
            if (token.progressBar != null) {
                token.progressBar.setValue(completedCount);
                token.progressBar.setString(completedCount + "/" + token.totalSubsamples);
            }
        });
    }

    @Override
    public void onSubsampleComplete(Object tokenObj, ArrayList<Double> stats, double Tn,
                                    double criticalValue, double pValue) {
        SubsampleToken token = (SubsampleToken) tokenObj;
        ArrayList<Double> copy = new ArrayList<>(stats);
        String rejectStr = Tn > criticalValue ? "REJECT" : "accept";
        boolean wasBottom = isScrolledToBottom();

        SwingUtilities.invokeLater(() -> {
            updateDensityPlot(token, copy, copy.size(), true, criticalValue, pValue);

            // Result label
            String resultText = String.format("    Tn=%.1f  CV(95%%)=%.1f  p=%.3f  %s",
                    Tn, criticalValue, pValue, rejectStr);
            JLabel resultLabel = new JLabel(resultText);
            resultLabel.setFont(new java.awt.Font("SansSerif", java.awt.Font.ITALIC, 10));
            resultLabel.setForeground(Tn > criticalValue ? new Color(180, 0, 0) : new Color(0, 130, 0));
            resultLabel.setAlignmentX(LEFT_ALIGNMENT);
            resultLabel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 16));

            // Fill progress bar and recolor
            if (token.progressBar != null) {
                int actualDraws = copy.size();
                token.progressBar.setValue(actualDraws);
                String doneLabel = actualDraws < token.totalSubsamples
                        ? "Done (B=" + actualDraws + "/" + token.totalSubsamples + ")"
                        : "Done";
                token.progressBar.setString(doneLabel);
                token.progressBar.setForeground(Tn > criticalValue
                        ? new Color(200, 60, 60) : new Color(50, 160, 50));
            }

            if (token.groupPanel != null) {
                token.groupPanel.add(resultLabel);
                token.groupPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, CHART_HEIGHT + 72));
                token.groupPanel.revalidate();
            }

            if (wasBottom) smartAutoScroll();
        });
    }

    // ---- KDE with Silverman reflection ----

    private void updateDensityPlot(SubsampleToken token, ArrayList<Double> data,
                                   int completedCount, boolean isFinal,
                                   double criticalValue, double pValue) {
        if (token.chart == null || data.isEmpty()) return;

        int n = data.size();
        double mean = data.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double variance = data.stream().mapToDouble(d -> (d - mean) * (d - mean)).sum() / Math.max(1, n - 1);
        double sd = Math.sqrt(variance);
        double bandwidth = 1.06 * sd * Math.pow(n, -0.2);
        if (bandwidth <= 0) bandwidth = 1.0;

        double minVal = 0; // bounded below by zero
        double maxVal = Collections.max(data) + 2 * bandwidth;
        maxVal = Math.max(maxVal, token.Tn + 2 * bandwidth);

        int numPoints = 150;
        XYSeries densitySeries = new XYSeries("Density");
        for (int i = 0; i <= numPoints; i++) {
            double x = minVal + (maxVal - minVal) * i / numPoints;
            double density = 0;
            for (double d : data) {
                double u = (x - d) / bandwidth;
                density += Math.exp(-0.5 * u * u);
                // Silverman reflection at zero
                double uR = (x + d) / bandwidth;
                density += Math.exp(-0.5 * uR * uR);
            }
            density /= (n * bandwidth * Math.sqrt(2 * Math.PI));
            densitySeries.add(x, density);
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(densitySeries);

        XYPlot plot = token.chart.getXYPlot();
        plot.setDataset(dataset);

        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
        Color lineColor;
        if (isFinal) {
            lineColor = token.Tn > criticalValue ? new Color(200, 60, 60) : new Color(50, 160, 50);
        } else {
            lineColor = new Color(50, 100, 180);
        }
        renderer.setSeriesPaint(0, lineColor);
        renderer.setSeriesStroke(0, new BasicStroke(isFinal ? 2.0f : 1.5f));
        plot.setRenderer(renderer);

        // Header label
        if (token.headerLabel != null) {
            String text;
            if (isFinal) {
                String result = token.Tn > criticalValue ? "REJECT" : "accept";
                text = "  " + token.prefix() + "Tree " + token.treeIndex + ", k=" + token.parameterIndex
                        + " \u2014 p=" + String.format("%.3f", pValue) + " " + result;
            } else {
                text = "  " + token.prefix() + "Tree " + token.treeIndex + ", k=" + token.parameterIndex
                        + " \u2014 " + completedCount + "/" + token.totalSubsamples;
            }
            token.headerLabel.setText(text);
        }

        // Markers
        plot.clearDomainMarkers();

        ValueMarker tnMarker = new ValueMarker(token.Tn);
        tnMarker.setPaint(Color.BLACK);
        tnMarker.setStroke(new BasicStroke(1.5f, BasicStroke.CAP_BUTT,
                BasicStroke.JOIN_MITER, 10.0f, new float[]{5.0f}, 0.0f));
        tnMarker.setLabel("Tn");
        tnMarker.setLabelFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 9));
        plot.addDomainMarker(tnMarker);

        if (isFinal && criticalValue > 0) {
            ValueMarker cvMarker = new ValueMarker(criticalValue);
            cvMarker.setPaint(Color.RED);
            cvMarker.setStroke(new BasicStroke(1.5f));
            cvMarker.setLabel("CV");
            cvMarker.setLabelFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 9));
            plot.addDomainMarker(cvMarker);
        } else if (data.size() >= 20) {
            ArrayList<Double> sorted = new ArrayList<>(data);
            Collections.sort(sorted);
            int idx95 = (int) Math.ceil(0.95 * sorted.size()) - 1;
            double running95 = sorted.get(Math.max(0, idx95));
            ValueMarker runningCV = new ValueMarker(running95);
            runningCV.setPaint(new Color(255, 100, 100, 180));
            runningCV.setStroke(new BasicStroke(1.0f));
            runningCV.setLabel("~95th");
            runningCV.setLabelFont(new java.awt.Font("SansSerif", java.awt.Font.PLAIN, 8));
            plot.addDomainMarker(runningCV);
        }
    }

    // ---- Proxy JTextArea that redirects to flow panel ----

    private static class FlowProxyTextArea extends JTextArea {
        private final SubsampleVisualizationPanel panel;

        FlowProxyTextArea(SubsampleVisualizationPanel panel) {
            this.panel = panel;
        }

        @Override
        public void append(String str) {
            if (str != null && !str.trim().isEmpty()) {
                panel.appendText(str);
            }
        }
    }
}
