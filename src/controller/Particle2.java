/*
package controller;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.Stack;
//import com.rits.cloning.Cloner;
import javafx.util.Pair;

public class Particle {

    private ArrayList<Double> position;
    private ArrayList<Double> pBest;
    static private Random rand = Main.rand;
    ArrayList<Double> velocity;
    private double fitness = 0;
    private double pBestFitness;
    //Cloner cloner = new Cloner();


    public ArrayList<Double> getVelocity() {
        return velocity;
    }

    public void setVelocity(ArrayList<Double> velocity) {
        this.velocity = velocity;
    }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public ArrayList<Double> getPosition() {
        return position;
    }

    public void setPosition(ArrayList<Double> position) {
        this.position = (ArrayList<Double>)position.clone();
        computeFitness();
    }

    public ArrayList<Double> getPBest() {
        return pBest;
    }

    public void setPBest(ArrayList<Double> pBest) {
        this.pBest = (ArrayList<Double>) pBest.clone();
        pBestFitness = fitness;
    }

    public Particle(int n){
        position = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double value = Main.round(Math.pow(2.5, rand.nextGaussian()), 1);
            if(value > 3.0)
                value = 3.0;
            position.add(value);
        }
        pBest = (ArrayList) position.clone();
        velocity  = new ArrayList<>();
        generateVelocity();
        computeFitness();
        pBestFitness = fitness;
    }

    public Particle (){
        position = new ArrayList<>();
        pBest = (ArrayList) position.clone();
        velocity  = new ArrayList<>();
        generateVelocity();
    }

    public ArrayList<Double[]> subtractPosition(ArrayList<Double> p1){
        ArrayList<Double[]> difference = new ArrayList<>();
        double incr = 0.05;
        for(int i = 0; i < p1.size(); i++){
            Double pVal = p1.get(i);
            if(!pVal.equals(position.get(i))){
//                double prevVelVal = velocity.get(i);
//                if(prevVelVal > pVal){
//                    velocity.set(i, prevVelVal - incr);
//                } else {
//                    velocity.set(i, prevVelVal + incr);
//                }
                Double[] value = new Double[2];
                value[0] = Double.valueOf(i);
                value[1] = p1.get(i);
                difference.add(value);
            }
        }
        return difference;
    }

    public void generateVelocity() {
        velocity.clear();
        for (int i = 0; i < position.size(); i++) {
            // Raidl and Bryant (2000)'s multiplicative scheme of bias based
            // normal distribution, with biasing strength set at 1.5 as they
            // recommended
            double bias = Main.round(Math.pow(2.5, rand.nextGaussian()), 1);

            // Clamp bias at 3.0
            if (bias > 3.0) {
                bias = 3.0;
            }
            velocity.add(bias);
        }
    }

    */
/**
     * Method which decodes the current position into its corresponding path
     * @return the cost of the constructed path
     *//*

    private double computeFitness(){
        */
/** First, we build the path according to position's biases **//*

        path = decodePath(position);
        */
/** Then, we compute its cost **//*

        fitness = getPathCost(path);
        return fitness;
    }

    public double getPathCost(Stack<Integer> path){
        double pathCost = 0;
        Iterator<Integer> i = path.iterator();
        int from = i.next();
        while (i.hasNext()){
            //System.out.println("in while");
            int to = i.next();
            //ArrayList<Pair<Double, Double>> costDistribution = Main.GRAPH.get(from).get(to);
            pathCost += Main.eCostGraph[from][to];
            from = to;
        }
        //System.out.println("out of while");
        return pathCost;
    }

    public static double getDistCost(ArrayList<Pair<Double, Double>> dist){
        double distCost = 0;
        //double prob = rand.nextDouble();
        for (int j = 0; j < dist.size(); j++) {
//            if(prob <= dist.get(j).getValue())
//                distCost = dist.get(j).getKey();
            //System.out.println(dist);
            distCost += dist.get(j).getKey() * dist.get(j).getValue();
        }
        //System.out.println("distCost: " + distCost);
        return distCost;
    }

    public static double getDistVariance(ArrayList<Pair<Double, Double>> dist){
        double expectedSum = 0;
        for (int j = 0; j < dist.size(); j++) {
            expectedSum += dist.get(j).getKey() * dist.get(j).getValue();
        }
        double expectedSquaredSum = 0;
        for (int j = 0; j < dist.size(); j++) {
            expectedSquaredSum += Math.pow(dist.get(j).getKey(), 2) * dist.get(j).getValue();
        }
        return expectedSquaredSum - Math.pow(expectedSum, 2);
    }


    public Stack<Integer> path;

    public Stack<Integer> decodePath(ArrayList<Double> biases){
        int start = Main.STARTNODE;
        int end = Main.ENDNODE;
        ArrayList<ArrayList<ArrayList<Pair<Double, Double>>>> graph = Main.GRAPH;
        path = new Stack<>();
        int currNode = start;
        path.push(currNode);
        ArrayList<Integer> invalid = new ArrayList<>();
//        for (int i = 0; i < Main.GRAPHSIZE; i++) {
//            graph[i][currNode] = 0;
//        }
        while(true){

            //System.out.println("inwhile: decodePath");
            double minBiasedCost = Double.MAX_VALUE;
            int nextNode = currNode;

            // Select node with least biased cost from adjacent nodes
            for (int i = 0; i < Main.GRAPHSIZE; i++) {
                if (invalid.contains(i) || path.contains(i) || Main.eCostGraph[currNode][i] == 0)
                    continue;
                // Bias edge cost by weight values of both nodes adjacent to it
                double biasedCost = biases.get(currNode) * biases.get(i)
                        * Main.eCostGraph[currNode][i];
                // Get edge with minimum biased cost
                if(biasedCost < minBiasedCost){
                    nextNode = i;
                    minBiasedCost = biasedCost;
                }
            }

            // Check if destination node has been reached
            if (nextNode == end){
                path.push(nextNode);
                break;
            }

            // Check for pedantic (dead end) node and backtrack if true
            if(minBiasedCost == Double.MAX_VALUE){
                int prevNode = path.pop();
                invalid.add(prevNode);
//                graph[prevNode][currNode] = 0;
//                graph[currNode][prevNode] = 0;
                currNode = path.peek();
                continue;
            }

            // Clear out the costs at nextNode's column of the adj matrix to
            // prevent duplication of nodes in path
//            for (int i = 0; i < graph.length; i++) {
//                graph[i][nextNode] = 0;
//            }
//            graph[nextNode][currNode] = 0;

            currNode = nextNode;
            path.push(currNode);
        }
        return path;
    }
    @Override
    public String toString(){
        return position.toString();
    }

    */
/**
     * @return the pBestFitness
     *//*

    public double getPBestFitness() {
        return pBestFitness;
    }

    */
/**
     * @param pBestFitness the pBestFitness to set
     *//*

    public void setPBestFitness(double pBestFitness) {
        this.pBestFitness = pBestFitness;
    }

    public void updatePBest(){
        pBest = (ArrayList) position.clone();
        pBestFitness = fitness;
    }

    */
/** @Override
    public Particle clone(){
    Particle p = new Particle(Main.GRAPHSIZE);
    p.setFitness(this.fitness);

    }
     **//*

}*/
