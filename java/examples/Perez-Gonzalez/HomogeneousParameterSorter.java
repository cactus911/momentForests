/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package examples.PerezGonzalez;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author stephen.p.ryan
 */
public class HomogeneousParameterSorter {

    public HomogeneousParameterSorter() {
    }
       
    public void sort(ArrayList<Integer> hpl, ArrayList<Double> hplStartingValues) {
        ArrayList<Duple> dupleList = new ArrayList<>();
        for(int i=0;i<hpl.size();i++) {
            dupleList.add(new Duple(hpl.get(i), hplStartingValues.get(i)));
        }
        Collections.sort(dupleList);
        hpl.clear();
        hplStartingValues.clear();
        for(Duple d : dupleList) {
            hpl.add(d.getIndex());
            hplStartingValues.add(d.getParameterValue());
        }
    }

    private class Duple implements Comparable<Duple> {

        int index;
        double parameterValue;

        public Duple(int index, double value) {
            this.index = index;
            this.parameterValue = value;
        }

        public int getIndex() {
            return index;
        }

        public double getParameterValue() {
            return parameterValue;
        }

        @Override
        public int compareTo(Duple o) {
            if (o.getIndex() > getIndex()) {
                return -1;
            }
            if (o.getIndex() < getIndex()) {
                return 1;
            }
            return 0;
        }

    }

}
