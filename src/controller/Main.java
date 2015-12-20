/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package controller;

import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;
import javafx.util.Pair;

import java.io.*;
import java.util.*;

//import com.rits.cloning.*;


/**
 *
 * @author Peter
 */

class Data implements Comparable<Data> {
    public final int index;
    public final double priority;


    public Data(int index, double priority) {
        this.index = index;
        this.priority = priority;
    }


    @Override
    public int compareTo(Data other) {
        return Double.valueOf(priority).compareTo(other.priority);
    }

    public boolean equals(Data other) {
        return priority == other.priority;
    }
}


public class Main {


    public static final int DP = 4;

    public static double c1 = 1.48;
    public static double c2 = 0.5;
    public static final int SEED = 12;
    public static Random rand = new Random();
    public static ArrayList<Particle> swarm;
    public static Particle dummyParticle;
    public static List<Particle> externalArchive;
    public static ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>> GRAPH;

    public static int GRAPHSIZE;
    public static int EDGE_NO;
    public static int STARTNODE;
    public static int ENDNODE;
    public static final int NO_OF_RUNS = 10;
    public static final int NO_OF_ITERATIONS = 50;
    public static final int SWARM_SIZE = 25;
    public static final double EPSILON = 0.01;
    public static final double epsPlus = EPSILON + 1;
    public static double minExpected;
    public static double minVariance;
    public static int k = 0;
    public static final double MUTATION_RATE = 0.2;
    public static ArrayList<String> improvedMinCosts = new ArrayList<>();
    public static ArrayList<String> improvedMinVars = new ArrayList<>();


    public static double eCostGraph[][];// = new double[GRAPHSIZE][GRAPHSIZE];
    public static double varianceGraph[][];// = new double[GRAPHSIZE][GRAPHSIZE];

    public static void main(String[] args) {
        dummyParticle = new Particle();
        long startTime = System.currentTimeMillis();
        psoSPP();
        long endTime = System.currentTimeMillis();
        double totalTime = (endTime-startTime)/1000.0;
        double averageTime = totalTime/50.0;
        System.out.println("\nTotal time: " + totalTime);
        System.out.println("Average time: " + averageTime);

    }

    public static void psoSPP(){
        try {
            PrintWriter pw = new PrintWriter(new File("result_1.txt"));

            double r1, r2;
            //GRAPH = randGraph();

            for (k = 0; k < NO_OF_RUNS; k++) {
                improvedMinCosts.add("New run");
                improvedMinVars.add("New run");
                System.out.println("Run " + (k+1) + ":");
                pw.println("Run " + (k+1) + ":");
                init();
                for (int i = 0; i < NO_OF_ITERATIONS; i++) {
                    for (Particle p : swarm) {
                        r1 = round(rand.nextDouble(), DP);
                        r2 = round(rand.nextDouble(), DP);

                        /** Selection of Leader **/
                        Particle leader = selectLeader(p);

                        // Get differences between X and pBest and X and gBest
                        ArrayList<Double[]> gDiff = p.subtractPosition(leader.getPosition());
                        ArrayList<Double[]> pDiff = p.subtractPosition(p.getPBest());

                        // Get magnitude of each difference
                        int pDiffMagnitude = (int) round(c1 * r1 * pDiff.size(), 0);
                        int gDiffMagnitude = (int) round(c2 * r2 * gDiff.size(), 0);

                        // Generate new velocity
                        Double[] newPosition = new Double[GRAPHSIZE];
                        Arrays.fill(newPosition, -1.0);
                        for (int j = 0; j < gDiffMagnitude; j++) {
                            int index = rand.nextInt(gDiff.size());
                            Double[] diffVal = gDiff.get(index);
                            newPosition[diffVal[0].intValue()] = diffVal[1];
                            gDiff.remove(index);
                            if(gDiff.isEmpty())
                                break;
                        }
                        for (int j = 0; j < pDiffMagnitude; j++) {
                            int index = rand.nextInt(pDiff.size());
                            Double[] diffVal = pDiff.get(index);
                            int pDiffIndex = diffVal[0].intValue();
                            if (newPosition[pDiffIndex] == -1.0) {
                                newPosition[pDiffIndex] = diffVal[1];
                            }
                            pDiff.remove(index);
                            if(pDiff.isEmpty())
                                break;
                        }
                        ArrayList<Integer> emptyIndices = new ArrayList<>();
                        for (int j = 0; j < newPosition.length; j++) {
                            if (newPosition[j] < 0) {
                                emptyIndices.add(j);
                            }
                        }

                        int prevVel = rand.nextInt(GRAPHSIZE) + 1;
                        int velCount = 0;
                        for (Integer j : emptyIndices) {
                            newPosition[j] = p.getVelocity().get(j);
                            velCount++;
                            if(emptyIndices.size() == GRAPHSIZE && velCount >= prevVel)
                                break;
                        }
                        for (int j = 0; j < newPosition.length; j++){
                            if(newPosition[j] < 0){
                                newPosition[j] = p.getPosition().get(j);
                            }
                        }


                        ArrayList<Double> newPositionList = new ArrayList<>(Arrays.asList(newPosition));

                        /** Mutation **/
                        if (rand.nextDouble() < MUTATION_RATE) {
                            int randIndex1 = rand.nextInt(newPositionList.size());
                            int randIndex2 = rand.nextInt(newPositionList.size());
                            double toMutate1 = newPositionList.get(randIndex1);
                            double toMutate2 = newPositionList.get(randIndex2);
                            if (toMutate1 <= 1.5) {
                                newPositionList.set(randIndex1, 2.9);
                                newPositionList.set(randIndex2, 0.1);
                            }
                            else {
                                newPositionList.set(randIndex1, 0.1);
                                newPositionList.set(randIndex2, 2.9);
                            }
                        }

                        /** Update of pBest **/
                        double prevVariance = Particle.getPathVar(Particle.decodePath(p.getPBest()));
                        double prevExpectedCost = Particle.getPathCost(Particle.decodePath(p.getPBest()));
                        p.setPosition(newPositionList);
                        p.generateVelocity();

                        if (p.getExpectedCost()/epsPlus <= prevExpectedCost || p.getVariance()/epsPlus <= prevVariance) {
                            p.updatePBest();
                            if(!isInExternalArchive(p))
                                externalArchive.add(new Particle(p));
                        }
                    }
                    updateExternalArchive();
                }
                double dijkstra = dijkstra(eCostGraph, STARTNODE, ENDNODE);
                double dijkstra2 = dijkstra(varianceGraph, STARTNODE, ENDNODE);

                System.out.println("Dijkstra cost: " + dijkstra);
                pw.println("Dijkstra cost: " + dijkstra);
                System.out.println("Dijkstra variance: " + dijkstra2);
                pw.println("Dijkstra variance: " + dijkstra2);
                System.out.println("Result: ");
                pw.println("Result: ");
                for (Particle p: externalArchive){
                    System.out.println(p.getPath() + ": Expected Cost = " + p.getExpectedCost() + " Variance = " + p.getVariance());
                    pw.println(p.getPath() + ": Expected Cost = " + p.getExpectedCost() + " Variance = " + p.getVariance());
                }
                System.out.println();
                pw.println();
            }

            System.out.println("\nImproved Min Costs");
            for (String s: improvedMinCosts){
                System.out.println(s);
            }

            System.out.println("\nImproved Min Variances");
            for (String s: improvedMinVars){
                System.out.println(s);
            }

            System.out.println();
            System.out.println("Mininmum Var = " + minVariance);
            System.out.println("Mininmum eCost = " + minExpected);
            pw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    private static boolean isInExternalArchive(Particle p) {
        for (Particle p1: externalArchive){
            if (p.getExpectedCost() == p1.getExpectedCost() && p.getVariance() == p1.getVariance())
                return true;
        }
        return false;
    }

    private static Particle selectLeader(Particle p) {

        double minDiff = Double.MAX_VALUE;
        int leaderIndex = 0;
        for (int i = 0; i < externalArchive.size(); i++) {
            Particle l = externalArchive.get(i);
            double lDiff = Math.abs(p.getSigmaValue() - l.getSigmaValue());
            if (lDiff < minDiff){
                minDiff = lDiff;
                leaderIndex = i;
            }
        }
        return externalArchive.get(leaderIndex);
    }

    public static void printVelocity(ArrayList<Double[]> v) {
        System.out.print("[");
        for (int i = 0; i < v.size() - 1; i++) {
            System.out.print(v.get(i)[0] + "(" + v.get(i)[1] + "), ");
        }
        System.out.println(v.get(v.size() - 1)[0] + "(" + v.get(v.size() - 1)[1] + ")]");
    }

    private static void init() {
        //GRAPH = randGraph();
        //System.out.println(GRAPH);
        readDistGraph();
        System.out.println(Arrays.deepToString(eCostGraph));
        minExpected = Double.MAX_VALUE;
        minVariance = Double.MAX_VALUE;
        swarm = new ArrayList<>();
        for (int i = 0; i < SWARM_SIZE; i++) {
            Particle p = new Particle(GRAPHSIZE);
            swarm.add(p);
        }
        System.out.println(swarm);
        computeNonDominated();
        for (Particle p: swarm)
            p.updateSigma();
        for (Particle p: externalArchive)
            p.updateSigma();
    }



    static double minFitness;

    private static void computeNonDominated(){
        externalArchive = new ArrayList<>();
        outerLoop: for (int i = 0; i < swarm.size(); i++) {
            Particle p1 = swarm.get(i);

            /** Update min objective values for sigma value calculation **/
            if (p1.getExpectedCost() < minExpected) {
                minExpected = p1.getExpectedCost();
                improvedMinCosts.add("iExpected Cost: " + p1.getExpectedCost() + " Variance: " + p1.getVariance());
            }
            if (p1.getVariance() < minVariance) {
                minVariance = p1.getVariance();
                improvedMinVars.add("iExpected Cost: " + p1.getExpectedCost() + " Variance: " + p1.getVariance());
            }


            for (int j = 0; j < swarm.size(); j++) {
                Particle p2 = swarm.get(j);
                if(p1.equals(p2))
                    continue;
                /*System.out.println("p1: " + p1);
                System.out.println(p1.getExpectedCost() + ", " + p1.getVariance());
                System.out.println("p2: " + p2);
                System.out.println(p2.getExpectedCost() + ", " + p2.getVariance());*/
                /** Checking for epsilon-non-dominance **/
                //System.out.println("Checked");
                if(p2.getExpectedCost() < p1.getExpectedCost() && p2.getVariance() < p1.getVariance()) {
                   continue outerLoop;
                }
            }
            /** Add epsilon-non-dominated candidate to external archive **/
            //System.out.println("Added");
            if(!isInExternalArchive(p1) || externalArchive.isEmpty()) {
                externalArchive.add(new Particle(p1));
                //ystem.out.println("Added " + p1 + " " + p1.getExpectedCost() + " " + p1.getVariance());
            }
        }
    }

    private static void updateExternalArchive(){
        outerLoop: for (int i = 0; i < externalArchive.size(); i++) {
            Particle p1 = externalArchive.get(i);

            /** Update min objective values for sigma value calculation **/
            if (p1.getExpectedCost() < minExpected) {
                minExpected = p1.getExpectedCost();
                improvedMinCosts.add("vExpected Cost: " + p1.getExpectedCost() + " Variance: " + p1.getVariance());
            }
            if (p1.getVariance() < minVariance) {
                minVariance = p1.getVariance();
                improvedMinVars.add("vExpected Cost: " + p1.getExpectedCost() + " Variance: " + p1.getVariance());
            }


            for (int j = 0; j < externalArchive.size(); j++) {
                Particle p2 = externalArchive.get(j);
                if(p1.equals(p2))
                    continue;
                /** Checking for and removing epsilon-dominated candidates **/
                if(p2.getExpectedCost()/epsPlus < p1.getExpectedCost() && p2.getVariance()/epsPlus < p1.getVariance()){
//                    System.out.println("Removing: " + p1 + " " + p1.getExpectedCost() + " " + p1.getVariance());
//                    System.out.println("Dominated by: " + p2 + " " + p2.getExpectedCost() + " " + p2.getVariance());
                    externalArchive.remove(p1);
                    i = 0;
                    continue outerLoop;
                } else if(p1.getExpectedCost()/epsPlus < p2.getExpectedCost() && p1.getVariance()/epsPlus < p2.getVariance()){
                    externalArchive.remove(p2);
                }
            }
        }
    }

    public static double round(double d, int numbersAfterDecimalPoint) {
        double n = Math.pow(10, numbersAfterDecimalPoint);
        double d2 = d * n;
        long lon = (long) d2;
        lon = ((long) (d2 + 0.5) > lon) ? lon + 1 : lon;
        return (lon) / n;
    }

    public static ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>> randGraph() {
        Random rand = new Random();
        ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>> g = new ArrayList<>(GRAPHSIZE);
        for (int i = 0; i < GRAPHSIZE ; i++) {
            ArrayList<ArrayList<Pair<Double, Double>>> row = new ArrayList<>();
            for (int j = 0; j < GRAPHSIZE; j++) {
                ArrayList<Pair<Double, Double>> pairList = new ArrayList<>();
                pairList.add(new Pair(-1.0, 1.0));
                row.add(pairList);
            }
            g.add(row);
        }
        int counter = 0;
        for (int i = 0; i < GRAPHSIZE; i++) {
            int j = rand.nextInt(GRAPHSIZE);
            if(i == j)
                continue;
            int distSize = 1 + rand.nextInt(5);
            ArrayList<Double> distribution = generateDistribution(distSize, 1000);
            ArrayList<Pair<Double, Double>> costProb = new ArrayList<>();
            for (int k = 0; k < distribution.size(); k++) {
                costProb.add(new Pair<>(distribution.get(k), 1.0 + rand.nextInt(1000)));
            }
            g.get(i).set(j, costProb);
            g.get(j).set(i, costProb);
            double expectedCost = Particle.getDistCost(costProb);
            double variance = Particle.getDistVariance(costProb);
            eCostGraph[i][j] = expectedCost;
            eCostGraph[j][i] = expectedCost;
            varianceGraph[i][j] = variance;
            varianceGraph[j][i] = variance;

            counter++;
            if (counter == EDGE_NO) {
                break;
            }
        }
        for (;;) {
            int i = rand.nextInt(GRAPHSIZE);
            int j = rand.nextInt(GRAPHSIZE);
            if (i == j) {
                continue;
            }
            if (g.get(i).get(j).get(0).getKey() < 0.0) {
                int distSize = 1 + rand.nextInt(5);
                ArrayList<Double> distribution = generateDistribution(distSize, 1000);
                ArrayList<Pair<Double, Double>> costProb = new ArrayList<>();
                for (int k = 0; k < distribution.size(); k++){
                    costProb.add(new Pair<>(distribution.get(k), 10.0 + rand.nextInt(1000)));
                }
                g.get(i).set(j, costProb);
                g.get(j).set(i, costProb);
                double expectedCost = Particle.getDistCost(costProb);
                double variance = Particle.getDistVariance(costProb);
                eCostGraph[i][j] = expectedCost;
                eCostGraph[j][i] = expectedCost;
                varianceGraph[i][j] = variance;
                varianceGraph[j][i] = variance;
                counter++;
            }
            if (counter == EDGE_NO) {
                break;
            }
        }
        return g;
    }

    @SuppressWarnings("Duplicates")
    public static ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>>  loadGraph(){
        ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>> g = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File("small_data_1.txt")));
            String s = br.readLine();
            System.out.println(s);
            GRAPHSIZE = Integer.parseInt(s);
            EDGE_NO = Integer.parseInt(br.readLine());
            String[] sHolder = br.readLine().split("\\s");
            STARTNODE = Integer.parseInt(sHolder[0]);
            ENDNODE = Integer.parseInt(sHolder[1]);
            g = new ArrayList<>(GRAPHSIZE);
            for (int i = 0; i < GRAPHSIZE; i++) {
                ArrayList<ArrayList<Pair<Double, Double>>> row = new ArrayList<>();
                for (int j = 0; j < GRAPHSIZE; j++) {
                    ArrayList<Pair<Double, Double>> pairList = new ArrayList<>();
                    pairList.add(new Pair<>(-1.0, 1.0));
                    row.add(pairList);
                }
                g.add(row);
            }
            while (br.ready()){
                sHolder = br.readLine().split("\\t");
                //System.out.println(Arrays.toString(sHolder));
                String[] iHolder = sHolder[0].split(",");
                int iIndex = Integer.parseInt(iHolder[0].substring(1));
                int jIndex = Integer.parseInt(iHolder[1].substring(0,iHolder[1].indexOf(")")));
                String[] cHolder = sHolder[1].split("\\s");
                String[] pHolder = sHolder[2].split("\\s");
                ArrayList<Pair<Double, Double>> pairList = new ArrayList<>();
                for (int i = 0; i < cHolder.length; i++) {
                    pairList.add(new Pair<>(Double.parseDouble(cHolder[i]), Double.parseDouble(pHolder[i])));
                }
//                System.out.println((iIndex-1) + ":" + (jIndex - 1));
//                System.out.println((jIndex-1) + ":" + (iIndex - 1));
//                System.out.println();
                g.get(iIndex-1).set(jIndex-1, pairList);
                g.get(jIndex-1).set(iIndex-1, pairList);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return g;
    }

    public static void readDistGraph(){
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File("rand_graph_dists.txt")));
            String s = br.readLine();
            System.out.println(s);
            String[] sHolder = s.split("\\s");
            GRAPHSIZE = Integer.parseInt(sHolder[0]);
            EDGE_NO = Integer.parseInt(sHolder[1]);
            eCostGraph = new double[GRAPHSIZE][GRAPHSIZE];
            varianceGraph = new double[GRAPHSIZE][GRAPHSIZE];

            sHolder = br.readLine().split("\\s");
            STARTNODE = Integer.parseInt(sHolder[0]);
            ENDNODE = Integer.parseInt(sHolder[1]);

            while (br.ready()) {
                s = br.readLine();
                s = s.replaceAll("[ \\n\\x0B\\f\\r]","");
                sHolder = s.split("\\t");
                //System.out.println(Arrays.toString(sHolder));
                String[] iHolder = sHolder[0].split(",");
                int iIndex = Integer.parseInt(iHolder[0].substring(1));
                int jIndex = Integer.parseInt(iHolder[1].substring(0, iHolder[1].indexOf(")")));

                double res[] = computeMeanVariance(sHolder[1]);
                eCostGraph[iIndex-1][jIndex-1] = res[0];
                varianceGraph[iIndex-1][jIndex-1] = res[1];
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] computeMeanVariance(String x){
        x = x.replaceAll("\\s+","");
        char dist = x.charAt(0);
        String s = x.substring(x.indexOf('(')+1, x.indexOf(')'));
        int a, b, c;
        double mean;
        double variance;
        switch(dist){
            case 'E':
                a = Integer.parseInt(s);
                mean = Math.pow(a, -1);
                variance = Math.pow(a, -1) * Math.log(2);
                return new double[] {mean, variance};
            case 'U':
                String[] sH = s.split(",");
                a = Integer.parseInt(sH[0]);
                b = Integer.parseInt(sH[1]);
                mean = 0.5 * (a + b);
                variance = (1.0/12) * Math.pow(b-a, 2);
                return new double[] {mean, variance};
            case 'N':
                sH = s.split(",");
                a = Integer.parseInt(sH[0]);
                b = Integer.parseInt(sH[1]);
                return new double[] {a, b};
            case 'T':
                sH = s.split(",");
                a = Integer.parseInt(sH[0]);
                b = Integer.parseInt(sH[1]);
                c = Integer.parseInt(sH[2]);
                mean = (a + b + c) / 3.0;
                variance = (Math.pow(a, 2) + Math.pow(b, 2) + Math.pow(c, 2) - (a * b)
                            - (a * c) - (b * c)) / 18.0;
                return new double[] {mean, variance};
            default:
                return null;

        }
    }

    /** Imported Dijkstra Algorithm with Priority Queue **/
    private static double dijkstra(double[][] G, int i, int j){
        //Get the number of vertices in G
        int n = G.length;
        //System.out.println(Arrays.deepToString(G));

		/* ... Your code here ... */
        double[] distance = new double[G.length];
        PriorityQueue<Data> PQ = new PriorityQueue<>();
        boolean[] inTree = new boolean[G.length];

        for (int index = 0; index < G.length; index++) {
            if (index == i) {
                distance[index] = 0;
            }
            else {
                distance[index] = Integer.MAX_VALUE;

                PQ.add(new Data(index, distance[index]));
                inTree[index] = true;
            }

        }

        for (int index = 0; index < G.length; index++) { // for each edge (v,z) do
            if (G[i][index] != 0) { // There is an edge
                if (distance[i] + G[i][index] < distance[index]) { // if D[v] + w((v,z)) < D[z] then
                    double oldIndex = distance[index];
                    distance[index] = distance[i] + G[i][index]; // D[z] ← D[v] + w((v,z))
                    Data t = new Data(index, oldIndex);
                    PQ.remove(new Data(index, oldIndex));
                    PQ.add(new Data(index, distance[index])); // update PQ wrt D[z]
                }
            }
        }


        while (PQ.peek() != null) { // If PQ isn't empty
            Data vertex = PQ.poll(); // RemoveMin
            for (int index = 0; index < G.length; index++) { // for each edge (u,z) with z ∈ PQ do
                if (G[vertex.index][index] != 0 && inTree[index]) { // z ∈ PQ
                    if (distance[vertex.index] + G[vertex.index][index] < distance[index]) { // if D[v] + w((v,z)) < D[z] then
                        double oldIndex = distance[index];
                        distance[index] = distance[vertex.index] + G[vertex.index][index]; // D[z] ← D[v] + w((v,z))
                        PQ.remove(new Data(index, oldIndex));
                        PQ.add(new Data(index, distance[index])); // update PQ wrt D[z]
                    }
                }

            }
        }
        if (distance[j] == Integer.MAX_VALUE || distance[j] < 0) {
            return -1;
        }
        else {
            return distance[j];
        }
    }
    private static ArrayList<Double> memeticSearch(ArrayList<Double> prevVal) {
        double min = minFitness;
        double gamma = 0.5;
        int prob = 0;
        ArrayList<Double> result = (ArrayList<Double>)prevVal.clone();
        ArrayList<Double> prospective = (ArrayList<Double>)prevVal.clone();
        ArrayList<Double> z = new ArrayList<>();
        for (int i = 0; i < GRAPHSIZE; i++) {
            z.add(rand.nextDouble() * 3);
        }
        //System.out.println("Previous: " + prevVal);

        for (int i = 0; i < 10; i++) {
            prevVal = (ArrayList) prospective.clone();
            int index = rand.nextInt(GRAPHSIZE);
            for (int j = 0; j < GRAPHSIZE; j++) {
                prob = rand.nextInt(2);
                //System.out.println(index);
                if(prob == 0)
                    gamma = gamma * -1;
                double newVal = prevVal.get(j) + gamma * z.get(j);
                if(newVal > 3.0)
                    newVal = 3.0;
                prospective.set(j, newVal);
            }
            //System.out.println(prospective);
            double fitness = dummyParticle.getPathCost(dummyParticle.decodePath(prospective));
            //System.out.println("Min: " + min + " MemFitness: " + fitness);
            if (fitness > min) {
                for (int j = 0; j < GRAPHSIZE; j++) {
                    z.set(j, rand.nextDouble() * 3);
                }
                gamma -= 0.01;
            }
            else if(fitness < min){
                min = fitness;
                System.out.println("Successfull memetic");
                result = (ArrayList<Double>) prospective.clone();
            }
        }
        return result;
    }

    public static ArrayList<Double> generateDistribution(int size, int max){
        Random rand = new Random();
        Set<Double> set = new TreeSet<>();
        set.add(0.0);
        while (set.size() < size) {
            set.add((1.0 + rand.nextInt(max)));
        }
        set.add((double) max);
        ArrayList<Double> output = new ArrayList<>();
        output.addAll(set);
        for (int i = 1; i < output.size(); i++) {
            double diff = output.get(i) - output.get(i-1);
            output.set(i-1, diff);
        }
        output.remove(output.size()-1);
        for (int i = 0; i < output.size(); i++) {
            double quotient = round(output.get(i)/max, 2);
            output.set(i, quotient);
        }
        return output;
    }
}