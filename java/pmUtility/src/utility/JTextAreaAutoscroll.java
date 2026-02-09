/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import javax.swing.JTextArea;

/**
 *
 * @author Stephen P. Ryan
 */
public class JTextAreaAutoscroll extends JTextArea {

    public JTextAreaAutoscroll() {
        super();
        setEditable(false);
    }

    public JTextAreaAutoscroll(String s) {
        super(s);
        setEditable(false);
    }
    
    @Override
    public void append(String str) {
        super.append(str);
        setCaretPosition(getDocument().getLength());
    }
}
