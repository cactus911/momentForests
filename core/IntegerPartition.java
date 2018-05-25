/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.util.ArrayList;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class IntegerPartition {

    ArrayList<Integer> left;
    ArrayList<Integer> right;

    public IntegerPartition(ArrayList<Integer> left, ArrayList<Integer> right) {
        this.left = left;
        this.right = right;
    }

    public ArrayList<Integer> getLeft() {
        return left;
    }

    public ArrayList<Integer> getRight() {
        return right;
    }

    @Override
    public String toString() {
        String s = "";
        s = s.concat("{ ");
        for (int i = 0; i < left.size(); i++) {
            s = s.concat(left.get(i) + " ");
        }
        s = s.concat("} { ");
        for (int i = 0; i < right.size(); i++) {
            s = s.concat(right.get(i) + " ");
        }
        s = s.concat("}");
        return s;
    }

}
