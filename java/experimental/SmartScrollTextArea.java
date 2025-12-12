/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package experimental;

import javax.swing.*;
import javax.swing.text.*;

public class SmartScrollTextArea extends JTextArea {
    private JScrollPane scrollPane;
    private DefaultCaret caret;
    
    public SmartScrollTextArea() {
        this(20, 50);
    }
    
    public SmartScrollTextArea(int rows, int columns) {
        super(rows, columns);
        setEditable(false);
        
        // Get the caret and set it to update policy ALWAYS_UPDATE
        caret = (DefaultCaret) getCaret();
        caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
    }
    
    /**
     * Call this after adding the JTextArea to a JScrollPane to enable smart scrolling
     */
    public void enableSmartScroll(JScrollPane scrollPane) {
        this.scrollPane = scrollPane;
        
        // Listen to document changes to check scroll position BEFORE new text
        getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            @Override
            public void insertUpdate(javax.swing.event.DocumentEvent e) {
                SwingUtilities.invokeLater(() -> updateScrollPolicy());
            }
            
            @Override
            public void removeUpdate(javax.swing.event.DocumentEvent e) {
                SwingUtilities.invokeLater(() -> updateScrollPolicy());
            }
            
            @Override
            public void changedUpdate(javax.swing.event.DocumentEvent e) {
                SwingUtilities.invokeLater(() -> updateScrollPolicy());
            }
        });
        
        // Also check on manual scrolling
        scrollPane.getVerticalScrollBar().addAdjustmentListener(e -> {
            if (!e.getValueIsAdjusting()) {
                updateScrollPolicy();
            }
        });
    }
    
    private void updateScrollPolicy() {
        if (scrollPane == null) return;
        
        JScrollBar scrollBar = scrollPane.getVerticalScrollBar();
        int extent = scrollBar.getModel().getExtent();
        int maximum = scrollBar.getModel().getMaximum();
        int value = scrollBar.getValue();
        
        // If we're at the bottom (within 1 pixel tolerance), enable auto-scroll
        if (value + extent >= maximum - 1) {
            caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
        } else {
            caret.setUpdatePolicy(DefaultCaret.NEVER_UPDATE);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Smart Scroll Demo");
            
            SmartScrollTextArea textArea = new SmartScrollTextArea();
            JScrollPane scrollPane = new JScrollPane(textArea);
            textArea.enableSmartScroll(scrollPane);
            
            frame.add(scrollPane);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.pack();
            frame.setVisible(true);
            
            // Simulate continuous text additions
            new Timer(100, e -> {
                textArea.append("New line at " + System.currentTimeMillis() + "\n");
            }).start();
        });
    }
}