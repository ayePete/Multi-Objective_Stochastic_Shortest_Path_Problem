/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package controller;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

/**
 *
 * @author Peter
 */
public class Test {

    public static Random rand = new Random();
    public static int limit = 30;
    public static int GRAPH_SIZE = 30;
    public static int EDGE_NO = 50;
    public static boolean inGraph[][] = new boolean[GRAPH_SIZE][GRAPH_SIZE];
    public static void main(String[] args){

                //System.out.println(Main.randGraph());
        /*try {
            System.out.println(Main.loadGraph());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }*/

        generateRandGraphDistributions();
    }

    public static void generateRandGraphDistributions(){
        try {
            int selectedDist = 0;
            int counter = 0;
            PrintWriter pw = new PrintWriter(new File("rand_graph_dists.txt"));
            for (int i = 0; i < GRAPH_SIZE; i++) {
                int j = rand.nextInt(GRAPH_SIZE);
                if(i == j) {
                    i--;
                    continue;
                }
                if(inGraph[i][j])
                    continue;

                selectedDist = rand.nextInt(3);
                pw.println(printGraph(i, j, selectedDist));
                inGraph[i][j] = true;

                counter++;
                if (counter == EDGE_NO) {
                    break;
                }
            }

            for (;;) {
                int i = rand.nextInt(GRAPH_SIZE);
                int j = rand.nextInt(GRAPH_SIZE);
                if (i == j) {
                    continue;
                }

                selectedDist = rand.nextInt(3);
                printGraph(i, j, selectedDist);
                counter++;
                if (counter == EDGE_NO) {
                    break;
                }
            }
            pw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static String generateRandNormDist(){
        int a = rand.nextInt(limit);
        int b = 1 + rand.nextInt(limit);
        String ret = "N(" +  a + ", " + b + ")";
        return ret;
    }

    public static String generateRandUniformDist(){
        int a = rand.nextInt(limit - 5);
        int b = rand.nextInt(limit);
        while (a >= b){
            b = rand.nextInt(limit);
        }
        String ret = "U(" +  a + ", " + b + ")";
        return ret;
    }

    public static String generateRandTriangularDist(){
        int a = limit;
        int b = rand.nextInt(limit);
        int c = rand.nextInt(limit);
        while (a >= b){
            a = 1 + rand.nextInt(limit);
        }

        while (c > b || a < c){
            c = rand.nextInt(limit);
        }
        String ret = "T(" +  a + ", " + b + ", " + c + ")";
        return ret;
    }

    public static String generateRandExponentialDist(){
        int a = 1 + rand.nextInt(limit);
        String ret = "E(" +  a + ")";
        return ret;
    }

    public static String printGraph(int i, int j, int selectedDist){
        i++;
        j++;
        String s = "";
        switch (selectedDist){
            case 0:
                s = "(" + i + ", " + j + ")\t" + generateRandNormDist();
                break;
            case 1:
                s = "(" + i + ", " + j + ")\t" + generateRandUniformDist();
                break;
            case 2:
                s = "(" + i + ", " + j + ")\t" + generateRandExponentialDist();
                break;
            case 3:
                s = "(" + i + ", " + j + ")\t" + generateRandTriangularDist();
                break;
            default:
                s = "Invalid";
                break;
        }
        System.out.println(s);
        return s;

    }
}
