package controller;

import javafx.util.Pair;

import java.util.*;


public class Particle {
	
    private ArrayList<Double> position;
    private ArrayList<Double> pBest;
    Random rand = Main.rand;
    ArrayList<Double> velocity;
    private double pBestFitness;
    private double sigmaValue;
    private double variance;
    private double expectedCost;
    public Stack<Integer> path;


    public Particle (){
        position = new ArrayList<>();
        pBest = (ArrayList) position.clone();
        velocity  = new ArrayList<>();
        generateVelocity();
    }
    
    public Particle(int n){
        position = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double value = Main.round(Math.pow(2.5, rand.nextGaussian()), 1);
            /*if(value > 3.0)
                value = 3.0;*/
            position.add(value);
        }
        pBest = (ArrayList) position.clone();
        velocity  = new ArrayList<>();
        path = decodePath(position);
        expectedCost = getPathCost();
        variance = getPathVar();
        updateSigma();
        generateVelocity();
    }

    public Particle(Particle p){
        position = (ArrayList<Double>) p.getPosition().clone();
        pBest = (ArrayList<Double>) p.getPBest().clone();
        velocity  = (ArrayList<Double>) p.getVelocity().clone() ;
        path = (Stack<Integer>) p.getPath().clone();
        expectedCost = p.getExpectedCost();
        variance = p.getVariance();
        sigmaValue = p.getSigmaValue();
    }
    
    public ArrayList<Double> getVelocity() {
        return velocity;
    }

    public void setVelocity(ArrayList<Double> velocity) {
        this.velocity = velocity;
    }

    public ArrayList<Double> getPosition() {
        return position;
    }

    public void setPosition(ArrayList<Double> position) {
        this.position = (ArrayList<Double>)position.clone();
        path = decodePath(position);
        expectedCost = getPathCost();
        variance = getPathVar();
        updateSigma();
    }

    public ArrayList<Double> getPBest() {
        return pBest;
    }

    public void setPBest(ArrayList<Double> pBest) {
        this.pBest = (ArrayList<Double>) pBest.clone();
    }
    
    public ArrayList<Double[]> subtractPosition(ArrayList<Double> p1){
        ArrayList<Double[]> difference = new ArrayList<>();
        double incr = 0.05;
        for(int i = 0; i < p1.size(); i++){
            Double pVal = p1.get(i);
            if(!pVal.equals(position.get(i))){
                Double[] value = new Double[2];
                value[0] = Double.valueOf(i);
                value[1] = p1.get(i);
                difference.add(value);
            }
        }
        return difference;
    }
    
    public final void generateVelocity() {
        velocity.clear();
        for (int i = 0; i < position.size(); i++) {
            // Raidl and Bryant (2000)'s multiplicative scheme of bias based 
            // normal distribution, with biasing strength set at 1.5 as they 
            // recommended
            double bias = Main.round(Math.pow(2.5, rand.nextGaussian()), 1);
            
            // Clamp bias at 3.0
            /*if (bias > 3.0) {
                bias = 3.0;
            }*/
            velocity.add(bias);
        }
    }
    
    /**
     * Method which decodes the current position into its corresponding path
     * @return the cost of the constructed path
     */
    public double getPathCost(){
        return getPathCost(path);
    }

    public double getPathVar(){
        return getPathVar(path);
    }

    public static double getPathCost(Stack<Integer> path){
        double pathCost = 0;
        Iterator<Integer> i = path.iterator();
        int from = i.next();
        while (i.hasNext()){
            int to = i.next();
            ArrayList<Pair<Double, Double>> costDistribution = Main.GRAPH.get(from).get(to);
            pathCost += getDistCost(costDistribution);
            from = to;
        }
        return pathCost;
    }

    public static double getPathVar(Stack<Integer> path){
        double pathCost = 0;
        Iterator<Integer> i = path.iterator();
        int from = i.next();
        while (i.hasNext()){
            int to = i.next();
            ArrayList<Pair<Double, Double>> costDistribution = Main.GRAPH.get(from).get(to);
            pathCost += getDistVariance(costDistribution);
            from = to;
        }
        return pathCost;
    }
    
    public static double getDistCost(ArrayList<Pair<Double, Double>> dist){
        double distCost = 0;
        for (int j = 0; j < dist.size(); j++) {
                distCost += dist.get(j).getKey() * dist.get(j).getValue();
            }
        return distCost;
    }

    public static double getDistVariance(ArrayList<Pair<Double, Double>> dist){
        double expectedSum = 0;
        for (int j = 0; j < dist.size(); j++) {
            expectedSum += dist.get(j).getKey() * dist.get(j).getValue();
        }
        double expectedSumSquared = Math.pow(expectedSum, 2);

        double expectedSquaredSum = 0;
        for (int j = 0; j < dist.size(); j++) {
            expectedSquaredSum += Math.pow(dist.get(j).getValue(), 2) * dist.get(j).getKey();
        }
        return expectedSquaredSum - expectedSumSquared;
    }

    
    public static Stack<Integer> decodePath(ArrayList<Double> biases){
        int start = Main.STARTNODE;
        int end = Main.ENDNODE;
        Stack<Integer> path = new Stack<>();
        int currNode = start;
        path.push(currNode);
        ArrayList<Integer> invalid = new ArrayList<>();
        double w1 = Math.abs(Math.sin(2 * Math.PI * Main.k / 200));
        double w2 = 1 - w1;
        while(true){

            //System.out.println("inwhile: decodePath");
            double minBiasedCost = Double.MAX_VALUE;
            int nextNode = currNode;

            // Select node with least biased cost from adjacent nodes
            for (int i = 0; i < Main.GRAPHSIZE; i++) {
                if (invalid.contains(i) || path.contains(i) || Main.GRAPH.get(currNode).get(i).get(0).getKey() < 0)
                    continue;
                // Bias edge cost by weight values of both nodes adjacent to it
//                double biasedCost = biases.get(currNode) * biases.get(i)
//                        * (w1*getDistCost(Main.GRAPH.get(currNode).get(i)) + w2 * getDistVariance(Main.GRAPH.get(currNode).get(i)));
                double biasedCost = biases.get(currNode) * biases.get(i)
                        * getDistCost(Main.GRAPH.get(currNode).get(i));
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

    public Stack<Integer> getPath() {
        return path;
    }

    /**
     * @return the pBestFitness
     */
    public double getPBestFitness() {
        return pBestFitness;
    }

    /**
     * @param pBestFitness the pBestFitness to set
     */
    public void setPBestFitness(double pBestFitness) {
        this.pBestFitness = pBestFitness;
    }

    public void updatePBest(){
        pBest = (ArrayList) position.clone();
    }

    public double getVariance() {
        return variance;
    }

    public double getSigmaValue() {
        return sigmaValue;
    }

    public double getExpectedCost() {
        return expectedCost;
    }

    public void updateSigma(){
        sigmaValue = (Math.pow((Main.minVariance * expectedCost), 2) - Math.pow((Main.minExpected * variance), 2))
                / (Math.pow((Main.minVariance * expectedCost), 2) + Math.pow((Main.minExpected * variance), 2));
    }

    public boolean equals (Object obj){
        if (!(obj instanceof Particle))
            return false;
        if (obj == this)
            return true;

        Particle rhs = (Particle) obj;
        return rhs.getPosition().equals(position);
    }

    /*public int hashCode(){
        return Objects.hashCode(position);
    }*/

    /** @Override
    public Particle clone(){
        Particle p = new Particle(Main.GRAPHSIZE);
        p.setFitness(this.fitness);

    }
    **/
}