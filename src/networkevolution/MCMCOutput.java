/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package networkevolution;

import java.math.BigDecimal;
import java.util.ArrayList;
import network.MetabolicNetwork;

/**
 *
 * @author mithani
 */
public class MCMCOutput {

    // variables to hold results for each iteration
    private ArrayList lstPaths;
    private Double[] insRates;
    private Double[] delRates;
    private BigDecimal[] pathProbabilities;
    private ArrayList lstNetworks;
    private Double[] dependenceProbabilities;
//    private ArrayList lstStartNetworks;
//    private ArrayList lstEndNetworks;
//    private ArrayList lstLeftPaths;
//    private ArrayList lstRightPaths;
//    private ArrayList leftPathProbabilities;
//    private ArrayList rightPathProbabilities;

    // other results
    private ArrayList lstVisitedNetworks;
    private int acceptCountPath;
    private int acceptCountParameter;
    private int acceptCountDependenceProbability;
//    private BigDecimal[][] likelihood;
//    private BigDecimal[] maxLikelihood;
//    private BigDecimal totalLikelihood;

    // Iteration Variables
    private MetabolicNetwork currentNetwork;
    private double[] currentRates;
//    private int[] currentPath;
//    private BigDecimal currentPathProbability;
//    private MetabolicNetwork currentStartNetwork;
//    private MetabolicNetwork currentEndNetwork;
//    private BigDecimal currentStartNetworkProbability;
    
    public MCMCOutput(int nIter) {
        this.lstPaths = new ArrayList(nIter + 1);
//        this.lstLeftPaths = new ArrayList(nIter + 1);
//        this.lstRightPaths = new ArrayList(nIter + 1);
//            this.lstStartNetworks = new ArrayList(nIter);
//            this.lstEndNetworks = new ArrayList(nIter);
        this.lstNetworks = new ArrayList(nIter + 1);
        this.insRates = new Double[nIter + 1];
        this.delRates = new Double[nIter + 1];
        this.pathProbabilities = new BigDecimal[nIter + 1];
        this.dependenceProbabilities = new Double[nIter + 1];
//        this.leftPathProbabilities = new ArrayList(nIter + 1);
//        this.rightPathProbabilities = new ArrayList(nIter + 1);

        this.lstVisitedNetworks = new ArrayList((int) (nIter * 0.4));
        this.acceptCountPath = 0;
        this.acceptCountParameter = 0;
    }

    // Getters / Setters   
    public Double[] getDelRates() {
        return delRates;
    }

    public void setDelRates(Double[] delRates) {
        this.delRates = delRates;
    }

    public Double[] getInsRates() {
        return insRates;
    }

    public void setInsRates(Double[] insRates) {
        this.insRates = insRates;
    }

    public Double[] getDependenceProbabilities() {
        return dependenceProbabilities;
    }

    public void setDependenceProbabilities(Double[] dependenceProbabilities) {
        this.dependenceProbabilities = dependenceProbabilities;
    }

//    public BigDecimal[][] getLikelihood() {
//        return likelihood;
//    }
//
//    public void setLikelihood(BigDecimal[][] likelihood) {
//        this.likelihood = likelihood;
//    }
//
//    public BigDecimal[] getMaxLikelihood() {
//        return maxLikelihood;
//    }
//
//    public void setMaxLikelihood(BigDecimal[] maxLikelihood) {
//        this.maxLikelihood = maxLikelihood;
//    }
//
//    public BigDecimal getTotalLikelihood() {
//        return totalLikelihood;
//    }
//
//    public void setTotalLikelihood(BigDecimal totalLikelihood) {
//        this.totalLikelihood = totalLikelihood;
//    }

//    public ArrayList getEndNetworks() {
//        return lstEndNetworks;
//    }
//
//    public void setEndNetworks(ArrayList lstEndNetworks) {
//        this.lstEndNetworks = lstEndNetworks;
//    }
    public ArrayList getPaths() {
        return lstPaths;
    }

    public void setPaths(ArrayList lstPaths) {
        this.lstPaths = lstPaths;
    }

//    public ArrayList getLeftPaths() {
//        return lstLeftPaths;
//    }
//
//    public void setLeftPaths(ArrayList lstLeftPaths) {
//        this.lstLeftPaths = lstLeftPaths;
//    }

    public ArrayList getNetworks() {
        return lstNetworks;
    }

    public void setNetworks(ArrayList lstNetworks) {
        this.lstNetworks = lstNetworks;
    }

//    public ArrayList getRightPaths() {
//        return lstRightPaths;
//    }
//
//    public void setRightPaths(ArrayList lstRightPaths) {
//        this.lstRightPaths = lstRightPaths;
//    }

    public MetabolicNetwork getCurrentNetwork() {
        return currentNetwork;
    }

    public void setCurrentNetwork(MetabolicNetwork currentNetwork) {
        this.currentNetwork = currentNetwork;
    }

//    public ArrayList getStartNetworks() {
//        return lstStartNetworks;
//    }
//
//    public void setStartNetworks(ArrayList lstStartNetworks) {
//        this.lstStartNetworks = lstStartNetworks;
//    }
    public ArrayList getVisitedNetworks() {
        return lstVisitedNetworks;
    }

    public void setVisitedNetworks(ArrayList lstVisitedNetworks) {
        this.lstVisitedNetworks = lstVisitedNetworks;
    }

    public BigDecimal[] getPathProbabilities() {
        return pathProbabilities;
    }

    public void setPathProbabilities(BigDecimal[] pathProbabilities) {
        this.pathProbabilities = pathProbabilities;
    }

//    public ArrayList getLeftPathProbabilities() {
//        return leftPathProbabilities;
//    }
//
//    public void setLeftPathProbabilities(ArrayList leftPathProbabilities) {
//        this.leftPathProbabilities = leftPathProbabilities;
//    }
//
//    public ArrayList getRightPathProbabilities() {
//        return rightPathProbabilities;
//    }
//
//    public void setRightPathProbabilities(ArrayList rightPathProbabilities) {
//        this.rightPathProbabilities = rightPathProbabilities;
//    }

    public int getAcceptCountPath() {
        return acceptCountPath;
    }

    public void setAcceptCountPath(int acceptCountPath) {
        this.acceptCountPath = acceptCountPath;
    }

    public int getAcceptCountParameter() {
        return acceptCountParameter;
    }

    public void setAcceptCountParameter(int acceptCountParameter) {
        this.acceptCountParameter = acceptCountParameter;
    }

    public int getAcceptCountDependenceProbability() {
        return acceptCountDependenceProbability;
    }

    public void setAcceptCountDependenceProbability(int acceptCountDependenceProbability) {
        this.acceptCountDependenceProbability = acceptCountDependenceProbability;
    }
    
//    public int[] getCurrentPath() {
//        return currentPath;
//    }
//
//    public void setCurrentPath(int[] currentPath) {
//        this.currentPath = currentPath;
//    }
//
//    public BigDecimal getCurrentPathProbability() {
//        return currentPathProbability;
//    }
//
//    public void setCurrentPathProbability(BigDecimal currentPathProbability) {
//        this.currentPathProbability = currentPathProbability;
//    }

    public double[] getCurrentRates() {
        return currentRates;
    }

    public void setCurrentRates(double[] currentRates) {
        this.currentRates = currentRates;
    }

//    public MetabolicNetwork getCurrentEndNetwork() {
//        return currentEndNetwork;
//    }
//
//    public void setCurrentEndNetwork(MetabolicNetwork currentEndNetwork) {
//        this.currentEndNetwork = currentEndNetwork;
//    }
//
//    public MetabolicNetwork getCurrentStartNetwork() {
//        return currentStartNetwork;
//    }
//
//    public void setCurrentStartNetwork(MetabolicNetwork currentStartNetwork) {
//        this.currentStartNetwork = currentStartNetwork;
//    }
//
//    public BigDecimal getCurrentStartNetworkProbability() {
//        return currentStartNetworkProbability;
//    }
//
//    public void setCurrentStartNetworkProbability(BigDecimal currentStartNetworkProbability) {
//        this.currentStartNetworkProbability = currentStartNetworkProbability;
//    }
//
}
