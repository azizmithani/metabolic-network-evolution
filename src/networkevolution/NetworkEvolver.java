/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package networkevolution;

import Jama.Matrix;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import network.DirectedReaction;
import network.MetabolicNetwork;
import systemobject.PhyloNode;
import utilities.Utilities;

/**
 *
 * @author mithani
 */
public class NetworkEvolver {

    static public final int EVOL_MODEL_INDEPENDENT_EDGE = 1;
    static public final int EVOL_MODEL_NEIGHBOUR_DEPENDENT = 2;
    static public final int EVOL_MODEL_HYBRID = 3;
    // default network size to be used for equilibrium probability approximation, if none specified
    static private int DEFAULT_SUB_NETWORK_SIZE = 4;

    // variable to store full networks indexed by network sequence
    HashMap lstNetworks = new HashMap();

    public double[] calculateHyperedgeProbabilities(MetabolicNetwork theNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int model, double insRate, double delRate, double dependenceProbability) {

        // calculate the edge probabilities according to the specified model
        double[] edgeRates = calculateHyperedgeChangeRates(theNetwork, coreEdges, prohibEdges, model, insRate, delRate, dependenceProbability);
        // get the sum of rates
        double totalRate = Utilities.sum(edgeRates);

        double[] edgeProbabilities = new double[edgeRates.length];
        if (totalRate == 0.0) {
            // each has equal probability
            for (int i = 0; i < edgeProbabilities.length; i++)
                edgeProbabilities[i] = 1.0 / edgeProbabilities.length;
        } else {
            // normalise the rates to get probabilities
            for (int i = 0; i < edgeProbabilities.length; i++)
                edgeProbabilities[i] = edgeRates[i] / totalRate;
        }
        return edgeProbabilities;
    }

    public double[] calculateHyperedgeChangeRates(MetabolicNetwork theNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int model, double insRate, double delRate, double dependenceProbability) {

        // calculate the rate matrix
        double[][] rateMatrix = getRateMatrix(insRate, delRate);
        // get the rates
        double insRateQ = rateMatrix[0][1];
        double delRateQ = rateMatrix[1][0];

        MetabolicNetwork refNetwork;
        if (model == EVOL_MODEL_NEIGHBOUR_DEPENDENT || model == EVOL_MODEL_HYBRID) {
            // create the reference network
            refNetwork = theNetwork.clone();
            refNetwork.markAllReactionsAsActive();
            refNetwork.markAllDirectedReactionsAsActive();
            refNetwork.setNodeEdgeIncidenceMatrix(refNetwork.createNodeEdgeIncidenceMatrix());
            refNetwork.setNeighbourhoodMatrix(refNetwork.createNeighbourhoodMatrix());

        } else
            refNetwork = null;

        // get the reaction sequence
        Byte[] rxnSequence = theNetwork.getReactionSequence();
        double[] edgeRates = new double[rxnSequence.length];
        int i = -1;
        Iterator it = theNetwork.getDirectedReactions().iterator();
        while (it.hasNext()) {
            // get the reaction object
            DirectedReaction reaction = (DirectedReaction) it.next();
            // increment the index
            i++;

            if (coreEdges.contains(reaction) || prohibEdges.contains(reaction))
                continue; // zero probability of changing for core & prohibited edges

            // set appropriate rate
            switch (model) {
                case EVOL_MODEL_INDEPENDENT_EDGE:
                    edgeRates[i] = hyperedgeChangeRate(theNetwork, null, null, i, insRateQ, delRateQ, 0.0);
                    break;
                case EVOL_MODEL_NEIGHBOUR_DEPENDENT:
                    edgeRates[i] = hyperedgeChangeRate(theNetwork, refNetwork, reaction, i, insRateQ, delRateQ, 1.0);
                    break;
                case EVOL_MODEL_HYBRID:
                    edgeRates[i] = hyperedgeChangeRate(theNetwork, refNetwork, reaction, i, insRateQ, delRateQ, dependenceProbability);
                    break;
            }
        } // end while it.hasNext()
        return edgeRates;

    }

    public double[][] calculateEdgeRateMatrix(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork, DirectedReaction reaction, int idx, int evolutionModel, double insRate, double delRate, double dependenceProbability) {
        // zero probability of changing for core & prohibited edges
        // get the change rate for the edge
        double edgeRate = calculateHyperedgeChangeRate(theNetwork, refNetwork, reaction, idx, evolutionModel, insRate, delRate, dependenceProbability);

        MetabolicNetwork complementNetwork = theNetwork.updateNetwork(reaction);
        double edgeRateComplement = calculateHyperedgeChangeRate(complementNetwork, refNetwork, reaction, idx, evolutionModel, insRate, delRate, dependenceProbability);

        // set correct rates
        double insRateQ;
        double delRateQ;
        if (reaction.isActive()) {
            //insRateQ = edgeRatesComplement[idx];
            insRateQ = edgeRateComplement;
            delRateQ = edgeRate;
        } else {
            insRateQ = edgeRate;
            //delRateQ = edgeRatesComplement[idx];
            delRateQ = edgeRateComplement;
        }
        // create rate matrix for this edge
        double[][] edgeRateMatrix = new double[][]{{-insRateQ, insRateQ}, {delRateQ, -delRateQ}};

        return edgeRateMatrix;
    }

    private double calculateHyperedgeChangeRate(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork, DirectedReaction reaction, int idx,
            int model, double insRateQ, double delRateQ, double dependenceProbability) {

        double edgeRate = 0.0;
        switch (model) {
            case EVOL_MODEL_INDEPENDENT_EDGE:
                edgeRate = hyperedgeChangeRate(theNetwork, null, null, idx, insRateQ, delRateQ, 0.0);
                break;
            case EVOL_MODEL_NEIGHBOUR_DEPENDENT:
                edgeRate = hyperedgeChangeRate(theNetwork, refNetwork, reaction, idx, insRateQ, delRateQ, 1.0);
                break;
            case EVOL_MODEL_HYBRID:
                edgeRate = hyperedgeChangeRate(theNetwork, refNetwork, reaction, idx, insRateQ, delRateQ, dependenceProbability);
                break;
            }
        return edgeRate;
    }

    private double hyperedgeChangeRate(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork,
            DirectedReaction reaction, int idx, double insRateQ, double delRateQ, double dependenceProbability) {

        // select appropriate rate
        double rate = theNetwork.getReactionSequence()[idx].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT) ? insRateQ : delRateQ;

        if (dependenceProbability > 0) { // dependenceProbabilty == 0 corresponds to independent edge model
            // calculate the neighbourhood weights 
            double neighWeight = calculateNeighbourhoodWeight(theNetwork, refNetwork, reaction, idx);

            rate *= (dependenceProbability * neighWeight + (1 - dependenceProbability));
        }
        // return the rate
        return rate;
    }

 private double hyperedgeChangeRateND(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork,
            DirectedReaction reaction, int idx, double insRateQ, double delRateQ) {

        // select appropriate rate
        double rate = theNetwork.getReactionSequence()[idx].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT) ? insRateQ : delRateQ;

        // calculate the neighbourhood weights
        double neighWeight = calculateNeighbourhoodWeight(theNetwork, refNetwork, reaction, idx);

        // get the rate
        rate *= neighWeight;
        // return the rate
        return rate;
    }

    private double[][] getRateMatrix(double insRate, double delRate) {

        double[][] rm = new double[][]{{-insRate, insRate}, {delRate, -delRate}};
//        double[][] one = new double[][]{{1.0, 1.0}, {1.0, 1.0}};
//
//        // get the matrices
//        Matrix A = new Matrix(rm);
//        Matrix B = new Matrix(one);
//
//        // solve the system (pi*Q = 0 or Q'*pi = 0)
//        Matrix X = A.transpose().plus(B).solve(B);
//
//        double[] pi = new double[]{X.get(0, 0), X.get(1, 1)};
//
//        for (int i = 0; i < rm.length; i++) {
//            double[] rmRow = rm[i];
//            for (int j = 0; j < rmRow.length; j++)
//                rmRow[j] *= pi[1 - i];
//        }
        return rm;

    }

    private double calculateNeighbourhoodWeight(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork, DirectedReaction reaction, int rxnIdx) {

        // get the list of active directed reations in the give network
        ArrayList activeReactions = theNetwork.getActiveDirectedReactions();

        // get the neighbours in the reference network
        HashMap neighbours = refNetwork.findReactionNeighbours(reaction, rxnIdx);

        ArrayList presentNeighbours = new ArrayList(neighbours.size());

        // find the neighbours that are present in the given network
        Iterator it = neighbours.values().iterator();
        while (it.hasNext()) {
            DirectedReaction neighbour = (DirectedReaction) it.next();

            if (activeReactions.contains(neighbour)) {
                // add to the list of present reactions
                presentNeighbours.add(neighbour);
            }
        }

        int presentCount = presentNeighbours.size();
        int totalEdges = theNetwork.getActiveDirectedReactions().size();
        // dont count the given reaction
        if (reaction.isActive()) {
            totalEdges--;
        }

        if (presentCount == 0 || totalEdges == 0) {
            presentCount = 1;
            totalEdges = theNetwork.getDirectedReactions().size() + 1;
        }

        double neighWeight = (double) presentCount / (double) totalEdges;


        return neighWeight;
    }

    private double[] calculateNeighbourhoodWeight_irreversible(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork, DirectedReaction reaction, int rxnIdx) {

        // get the list of active directed reations in the give network
        ArrayList activeReactions = theNetwork.getActiveDirectedReactions();

        // get the neighbours in the reference network
        HashMap neighbours = refNetwork.findReactionNeighbours(reaction, rxnIdx);

        ArrayList presentNeighbours = new ArrayList(neighbours.size());
        ArrayList absentNeighbours = new ArrayList(neighbours.size());

        // keep track of number of reversible reactions
        int reversibleCountPresent = 0;
        int reversibleCountAbsent = 0;
        // find the neighbours that are present in the given network
        Iterator it = neighbours.values().iterator();
        while (it.hasNext()) {
            DirectedReaction neighbour = (DirectedReaction) it.next();

            if (activeReactions.contains(neighbour)) {
                // add to the list of present reactions
                presentNeighbours.add(neighbour);

                if (neighbour.isReversible())
                    reversibleCountPresent++;
            } else {
                // add to the list of absent reactions
                absentNeighbours.add(neighbour);

                if (neighbour.isReversible())
                    reversibleCountAbsent++;

            }
        }

        int presentCount = presentNeighbours.size();
        int absentCount = absentNeighbours.size();
        int totalPresent = theNetwork.getActiveDirectedReactions().size();
        int totalAbsent = theNetwork.getDirectedReactions().size() - totalPresent;

        // dont count the given reaction
        if (reaction.isActive()) {
            totalPresent--;
        } else {
            totalAbsent--;
        }

        // handle the case resulting in 0 weight, i.e. when all neighbours are either present or absent (or no neighbours at all)
        if (presentCount == 0) {
            presentCount = 1;
            totalPresent = theNetwork.getDirectedReactions().size();
        }
        if (absentCount == 0) {
            absentCount = 1;
            totalAbsent = theNetwork.getDirectedReactions().size();
        }

        double[] neighWeight = new double[2];
        neighWeight[0] = (totalPresent == 0 ? 1 : (double) presentCount / (double) totalPresent);
        neighWeight[1] = (totalAbsent == 0 ? 1 : (double) absentCount / (double) totalAbsent);

        return neighWeight;
    }

    static public void summariseEvolution(MetabolicNetwork mtbNetwork, ArrayList lstNetworks) {

        Byte[] theNetwork = null;
        Iterator it = lstNetworks.iterator();
        // get the first network
        if (it.hasNext())
            theNetwork = (Byte[]) it.next();
        else
            return;

        // variables to hold summary data 
        int[] insFreq = new int[theNetwork.length]; // number of insertions for each edge
        int[] delFreq = new int[theNetwork.length]; // number of deletions for each edge
        int[] edgeFreq = new int[theNetwork.length]; // number of times an edge is present
        int[] edgeCount = new int[lstNetworks.size()]; // number of edges present at each iteration


        summariseEvolution(lstNetworks, insFreq, delFreq, edgeFreq, edgeCount);

        // print the summary
        System.out.println("Edge No." + "\t" + "ID" + "\t" + "Ins Freq" + "\t" + "Del Freq" + "\t" + "Edge Freq");
        int idx = 0;
        it = mtbNetwork.getDirectedReactions().iterator();
        while (it.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) it.next();
            System.out.print(idx);
            System.out.print("\t" + reaction.getID());
            System.out.print("\t" + Integer.toString(insFreq[idx]));
            System.out.print("\t" + Integer.toString(delFreq[idx]));
            System.out.print("\t" + Integer.toString(edgeFreq[idx]));

            System.out.println();
            idx++;
        }

        System.out.println("Edge Counts");
        for (int i = 0; i < edgeCount.length; i++) {
            System.out.print(i);
            System.out.print("\t" + Integer.toString(edgeCount[i]));

            System.out.println();
        }

    }

    static public void summariseEvolution(ArrayList lstNetworks,
            int[] insFreq, int[] delFreq, int[] edgeFreq, int[] edgeCount) {

        Byte[] theNetwork = null, prevNetwork = null;
        Iterator it = lstNetworks.iterator();
        // get the first network
        if (it.hasNext())
            prevNetwork = (Byte[]) it.next();
        else
            return;

        if (edgeCount != null)
            // number of edge present for first network
            edgeCount[0] = Utilities.sum(prevNetwork);

        int iter = 1;
        while (it.hasNext()) {
            // get the next network
            theNetwork = (Byte[]) it.next();

            for (int i = 0; i < theNetwork.length; i++) {
                if (edgeFreq != null) {
                    // edge freq
                    if (theNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                        edgeFreq[i]++;
                }
                // ins / del freq
                if (!theNetwork[i].equals(prevNetwork[i])) {
                    if (insFreq != null && theNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) // edge was inserted
                        insFreq[i]++;
                    else if (delFreq != null && theNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT)) // edge was deleted
                        delFreq[i]++;
                }

            }

            if (edgeCount != null)
                // number of edge present for this iteration
                edgeCount[iter++] = Utilities.sum(theNetwork);

            prevNetwork = theNetwork;
        } // while it.hasNext()

    }

    public ArrayList allPathsK(int K, MetabolicNetwork startNetwork, MetabolicNetwork endNetwork,
            ArrayList coreEdges, ArrayList prohibEdges) {

        // setup the reaction sequences if not already present
        if (startNetwork.getReactionSequence() == null)
            startNetwork.setupNetworkSequence();
        if (endNetwork.getReactionSequence() == null)
            endNetwork.setupNetworkSequence();

        // initialise the variables
        ArrayList lstPaths = new ArrayList();
        ArrayList partialPath = new ArrayList(K);
        // call the recursive function
        allPathsK(startNetwork, endNetwork, coreEdges, prohibEdges, lstPaths, partialPath, K, K);

        lstPaths.trimToSize();

        return lstPaths;
    }

    private void allPathsK(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, ArrayList lstPaths,
            ArrayList partialPath, int remainingSteps, int K) {

        if (remainingSteps == 0)
            return;
        else {
            // Get the edges that can be inserted into the given network (inactive edges - prohibited edges)
            ArrayList insertableEdges = startNetwork.getInactiveDirectedReactions();
            insertableEdges.removeAll(prohibEdges);

            Iterator it = insertableEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction edge = (DirectedReaction) it.next();
                // update the network
                MetabolicNetwork theNetwork = startNetwork.updateNetwork(edge);

                // update the partial path
                ArrayList updatedPath = (ArrayList) partialPath.clone();
                updatedPath.add(startNetwork.getDirectedReactions().indexOf(edge));
                updatedPath.trimToSize();

                // if we have reached the final network in exact K steps then add the path to the list
                if (remainingSteps == 1 && Arrays.equals(theNetwork.getReactionSequence(), endNetwork.getReactionSequence()))
                    lstPaths.add(updatedPath.toArray());
                else
                    // otherwise generate the recursive call
                    allPathsK(theNetwork, endNetwork, coreEdges, prohibEdges, lstPaths,
                            updatedPath, remainingSteps - 1, K);
            } // while it.hasNext()

            // Get the edges that can be deleted from the given network (active edges - core edges)
            ArrayList deletableEdges = startNetwork.getActiveDirectedReactions();
            deletableEdges.removeAll(coreEdges);

            it = deletableEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction edge = (DirectedReaction) it.next();
                // update the network
                MetabolicNetwork theNetwork = startNetwork.updateNetwork(edge);

                // update the partial path
                ArrayList updatedPath = (ArrayList) partialPath.clone();
                updatedPath.add(startNetwork.getDirectedReactions().indexOf(edge));
                updatedPath.trimToSize();

                // if we have reached the final network in exact K steps then add the path to the list
                if (remainingSteps == 1 && Arrays.equals(theNetwork.getReactionSequence(), endNetwork.getReactionSequence()))
                    lstPaths.add(updatedPath.toArray());
                else
                    // otherwise generate the recursive call
                    allPathsK(theNetwork, endNetwork, coreEdges, prohibEdges, lstPaths,
                            updatedPath, remainingSteps - 1, K);
            } // while it.hasNext()

        } // end if remaining steps = 0
    }

    static public HashMap pathLengthDistribution(ArrayList pathList, int nBurning, boolean printOutput) {

        int[] pathLengths = new int[pathList.size() - nBurning - 1]; // ignore the starting point + burning period

        // get the path length. Also find the maximum & minimum path length
        int maxLen = 0;
        int minLen = Integer.MAX_VALUE;
        int itIdx = 0;
        // process all paths 
        Iterator it = pathList.iterator();
        while (it.hasNext()) {
            int[] path = (int[]) it.next();
            if (itIdx++ <= nBurning)
                continue;

            pathLengths[itIdx - nBurning - 2] = path.length;
            // check if it's the max length so far
            if (path.length > maxLen)
                maxLen = path.length;
            else if (path.length < minLen) // check if it's the minlength so far
                minLen = path.length;

        }

        // calculate the distribution of path lengths
        int[] pathLength = new int[maxLen + 1]; // ignore the 0th index
        for (int i = 0; i < pathLengths.length; i++) {
            pathLength[pathLengths[i]]++;
        }

        if (printOutput)
            System.out.println("Path Length Distribution:");
        HashMap pathLengthDistribution = new HashMap(pathLength.length / 2); // ignore the 0th index
        for (int i = minLen; i < pathLength.length; i = i + 2) {
            Integer length = new Integer(i);
            Integer freq = pathLength[i];
            pathLengthDistribution.put(length, freq);

            if (printOutput)
                System.out.println(length.toString() + "\t" + freq.toString());
        }


        return pathLengthDistribution;
    }

    static public HashMap pathDistribution(ArrayList pathList, int nBurning, boolean printOutput) {


        int[] pathFreq = new int[pathList.size() - nBurning - 1]; // ignore the starting point + burning period

        ArrayList processedPaths = new ArrayList((int) (pathList.size() * 0.4));
        ArrayList distinctPaths = new ArrayList((int) (pathList.size() * 0.4));
        int pathIdx = 0;
        int itIdx = 0;
        // process all paths 
        Iterator it = pathList.iterator();
        while (it.hasNext()) {
            int[] path = (int[]) it.next();
            if (itIdx++ <= nBurning)
                continue;

            String strPath = Utilities.toString(path);
            int idx = processedPaths.indexOf(strPath);
            if (idx == -1) {
                // we haven't seen this path before .. add it
                processedPaths.add(strPath);
                distinctPaths.add(path);
                pathFreq[pathIdx++] = 1;
            } else // increment the frequency
                pathFreq[idx]++;
        }

        if (printOutput)
            System.out.println("Path Distribution:");

        HashMap pathDistribution = new HashMap(distinctPaths.size());
        // process all distinct paths
        int idx = 0;
        it = distinctPaths.iterator();
        while (it.hasNext()) {
            int[] path = (int[]) it.next();
            Integer freq = new Integer(pathFreq[idx++]);
            pathDistribution.put(path, freq);
            if (printOutput)
                System.out.println(Utilities.toString(path) + "\t" + freq.toString());
        }

        return pathDistribution;
    }

    static public HashMap networkDistribution(MetabolicNetwork startNetwork, ArrayList distinctPathList, boolean printOutput) {

        ArrayList lstEdges = startNetwork.getDirectedReactions();

        int[] networkFreq = new int[distinctPathList.size() * ((int[]) distinctPathList.get(0)).length];

        ArrayList processedNetworks = new ArrayList((int) (distinctPathList.size()));
        ArrayList distinctNetworks = new ArrayList((int) (distinctPathList.size()));
        int networkIdx = 0;
        // process all paths 
        Iterator it = distinctPathList.iterator();
        while (it.hasNext()) {
            int[] path = (int[]) it.next();

            MetabolicNetwork theNetwork = startNetwork.clone();
            for (int j = 0; j < path.length - 1; j++) { // dont need to process the last event. it will always lead to destination network
                // get the edge corresponding to ith event
                DirectedReaction reaction = (DirectedReaction) lstEdges.get(path[j] - 1);

                // Update the network
                theNetwork = theNetwork.updateNetwork(reaction);

                // If we have not seen this network before then add it
                // to the list
                String strNetwork = Utilities.toString(theNetwork.getReactionSequence());
                int idx = processedNetworks.indexOf(strNetwork);
                if (idx == -1) {
                    // we haven't seen this network before .. add it
                    processedNetworks.add(strNetwork);
                    distinctNetworks.add(theNetwork.getReactionSequence());
                    networkFreq[networkIdx++] = 1;
                } else // increment the frequency
                    networkFreq[idx]++;
            } // for each event int the path
        } // for each distinct path

        if (printOutput)
            System.out.println("Network Distribution:");

        HashMap networkDistribution = new HashMap(distinctNetworks.size());
        // process all distinct paths
        int idx = 0;
        it = distinctNetworks.iterator();
        while (it.hasNext()) {
            Byte[] network = (Byte[]) it.next();
            Integer freq = new Integer(networkFreq[idx++]);
            networkDistribution.put(network, freq);

            if (printOutput)
                System.out.println(Utilities.toString(network) + "\t" + freq.toString());
        }

        return networkDistribution;
    }

    static public HashMap rateDistribution(Double[] rates, int nBurning, boolean printOutput, String heading) {


        int[] rateFreq = new int[rates.length - nBurning - 1]; // ignore the starting point + burning period

        ArrayList distinctRates = new ArrayList((int) (rates.length * 0.4));
        int rateIdx = 0;
        // process all paths 
        for (int i = nBurning + 1; i < rates.length; i++) {
            Double rate = rates[i];
            int idx = distinctRates.indexOf(rate);
            if (idx == -1) {
                // we haven't seen this rate before .. add it
                distinctRates.add(rate);
                rateFreq[rateIdx++] = 1;
            } else // increment the frequency
                rateFreq[idx]++;
        }


        if (printOutput && heading != null)
            System.out.println(heading);

        HashMap rateDistribution = new HashMap(distinctRates.size());
        // process all distinct rates
        int idx = 0;
        Iterator it = distinctRates.iterator();
        while (it.hasNext()) {
            Double rate = (Double) it.next();
            Integer freq = new Integer(rateFreq[idx++]);
            rateDistribution.put(rate, freq);

            if (printOutput && rate != null && freq != null)
                System.out.println(rate.toString() + "\t" + freq.toString());
        }

        return rateDistribution;
    }

    public BigDecimal[][] generateRateMatrix(ArrayList networkList, MetabolicNetwork refNetwork, int evolutionModel, BigDecimal insRate, BigDecimal delRate, double dependenceProbability) {

        if (networkList == null)
            return null;

        if (networkList.isEmpty())
            return new BigDecimal[0][0];

        // calculate the rate matrix
        double[][] edgeRateMatrix = getRateMatrix(insRate.doubleValue(), delRate.doubleValue());
        // get the rates
        double insRateQ = edgeRateMatrix[0][1];
        double delRateQ = edgeRateMatrix[1][0];

        // and the network count (from the network list)
        int nwCount = networkList.size();
        BigDecimal[][] rateMatrix = new BigDecimal[nwCount][nwCount];

        // initialise the rate matrix
        for (int i = 0; i < rateMatrix.length; i++) {
            BigDecimal[] row = rateMatrix[i];
            for (int j = 0; j < row.length; j++) {
                row[j] = BigDecimal.ZERO;
            }
        }

        int idxI = -1;
        Iterator itX = networkList.iterator();
        while (itX.hasNext()) {
            // get the next network sequence
            Byte[] startNetworkSeq = (Byte[]) itX.next();
            idxI++; // increment the index

            // get the corresponding network
            MetabolicNetwork startNetwork = refNetwork.getNetwork(startNetworkSeq);
            int idxJ = idxI;
            List partialList = networkList.subList(idxI + 1, nwCount);
            Iterator itY = partialList.iterator();
            while (itY.hasNext()) {
                // get the next network sequence
                Byte[] endNetworkSeq = (Byte[]) itY.next();
                idxJ++; // increment the index

                // get the corresponding network
                MetabolicNetwork endNetwork = refNetwork.getNetwork(endNetworkSeq);

                ArrayList lstDifferences = MetabolicNetwork.findDifferences(startNetwork, endNetwork);
                if (lstDifferences.size() == 1) {// the networks differ by only one edge

//                    System.out.println(Integer.toString(idxI) + ", " + Integer.toString(idxJ));
                    DirectedReaction reaction = (DirectedReaction) lstDifferences.get(0);
                    rateMatrix[idxI][idxJ] = new BigDecimal(calculateHyperedgeChangeRate(startNetwork, refNetwork, reaction, startNetwork.getDirectedReactions().indexOf(reaction), evolutionModel, insRateQ, delRateQ, dependenceProbability));
                    rateMatrix[idxJ][idxI] = new BigDecimal(calculateHyperedgeChangeRate(endNetwork, refNetwork, reaction, startNetwork.getDirectedReactions().indexOf(reaction), evolutionModel, insRateQ, delRateQ, dependenceProbability));
//                    if (startNetworkSize < endNetworkSize) {
//                        // toNetwork can be reached by an insertion event
//                        rateMatrix[idxI][idxJ] = insRate;
//                        // and the reverse case
//                        rateMatrix[idxJ][idxI] = delRate;
//                    } else {
//                        // toNetwork can be reached by a deletion event
//                        rateMatrix[idxI][idxJ] = delRate;
//                        // and the reverse case
//                        rateMatrix[idxJ][idxI] = insRate;
//                    }
                } // networks differ by single edge
            } // end for j
        } // end for i

        // set the diagonal elements
        for (int i = 0; i < rateMatrix.length; i++) {
            BigDecimal totRate = Utilities.sum(rateMatrix[i]);
            rateMatrix[i][i] = totRate.negate();
        }

        return rateMatrix;
    }

// Modified to use hashing to avoid calls to MetabolicNetwork.getNetwork()
//    public double[][] generateRateMatrix(ArrayList networkList, MetabolicNetwork refNetwork, int evolutionModel, double insRate, double delRate) {
//
//        if (networkList == null)
//            return null;
//
//        if (networkList.isEmpty())
//            return new double[0][0];
//
//        // calculate the rate matrix
//        double[][] edgeRateMatrix = getRateMatrix(insRate, delRate);
//        // get the rates
//        double insRateQ = edgeRateMatrix[0][1];
//        double delRateQ = edgeRateMatrix[1][0];
//
//        // and the network count (from the network list)
//        int nwCount = networkList.size();
//        double[][] rateMatrix = new double[nwCount][nwCount];
//
//        // initialise the rate matrix
//        for (int i = 0; i < rateMatrix.length; i++) {
//            double[] row = rateMatrix[i];
//            for (int j = 0; j < row.length; j++) {
//                row[j] = 0.0;
//            }
//        }
//
//        int idxI = -1;
//        Iterator itX = networkList.iterator();
//        while (itX.hasNext()) {
//            // get the next network sequence
//            Byte[] startNetworkSeq = (Byte[]) itX.next();
//            idxI++; // increment the index
//
//            // get the corresponding network
//            MetabolicNetwork startNetwork = refNetwork.getNetwork(startNetworkSeq);
//            // get a local copy of directed reactions
//            ArrayList lstDirectedReactionsStart = startNetwork.getDirectedReactions();
//            Object[] arrDirectedReactionsStart = startNetwork.getDirectedReactions().toArray();
//
//            int idxJ = idxI;
//            List partialList = networkList.subList(idxI + 1, nwCount);
//            Iterator itY = partialList.iterator();
//            while (itY.hasNext()) {
//                // get the next network sequence
//                Byte[] endNetworkSeq = (Byte[]) itY.next();
//                idxJ++; // increment the index
//
//                // get the corresponding network
//                MetabolicNetwork endNetwork = refNetwork.getNetwork(endNetworkSeq);
//                // get a local copy of directed reactions
//                ArrayList lstDirectedReactionsEnd = endNetwork.getDirectedReactions();
//
//                ArrayList lstDifferences = MetabolicNetwork.findDifferences(startNetwork, endNetwork);
//                if (lstDifferences.size() == 1) {// the networks differ by only one edge
//                    int rxnIdxStart = lstDirectedReactionsStart.indexOf(lstDifferences.get(0));
//                    int rxnIdxEnd = lstDirectedReactionsEnd.indexOf(lstDifferences.get(0));
//                    rateMatrix[idxI][idxJ] = calculateHyperedgeChangeRate(startNetwork, refNetwork, (DirectedReaction) arrDirectedReactionsStart[rxnIdxStart], rxnIdxStart, evolutionModel, insRateQ, delRateQ);
////                    rateMatrix[idxI][idxJ] = calculateHyperedgeChangeRate(startNetwork, refNetwork, (DirectedReaction) lstDirectedReactionsStart.get(rxnIdxStart), rxnIdxStart, evolutionModel, insRateQ, delRateQ);
//                    rateMatrix[idxJ][idxI] = calculateHyperedgeChangeRate(endNetwork, refNetwork, (DirectedReaction) lstDirectedReactionsEnd.get(rxnIdxEnd), rxnIdxEnd, evolutionModel, insRateQ, delRateQ);
//                } // networks differ by single edge
//            } // end for j
//        } // end for i
//
//        // set the diagonal elements
//        for (int i = 0; i < rateMatrix.length; i++) {
//            double totRate = Utilities.sum(rateMatrix[i]);
//            rateMatrix[i][i] = -totRate;
//        }
//
//        return rateMatrix;
//    }
    public double[][] generateRateMatrix(ArrayList networkList, MetabolicNetwork refNetwork, int evolutionModel, double insRate, double delRate, double dependenceProbability) {

        String strNetwork;

        if (networkList == null)
            return null;

        if (networkList.isEmpty())
            return new double[0][0];

        // calculate the rate matrix
        double[][] edgeRateMatrix = getRateMatrix(insRate, delRate);
        // get the rates
        double insRateQ = edgeRateMatrix[0][1];
        double delRateQ = edgeRateMatrix[1][0];

        // and the network count (from the network list)
        int nwCount = networkList.size();
        double[][] rateMatrix = new double[nwCount][nwCount];

        // initialise the rate matrix
        for (int i = 0; i < rateMatrix.length; i++) {
            double[] row = rateMatrix[i];
            for (int j = 0; j < row.length; j++) {
                row[j] = 0.0;
            }
        }

        int idxI = -1;
        Iterator itX = networkList.iterator();
        while (itX.hasNext()) {
            // get the next network sequence
            Byte[] startNetworkSeq = (Byte[]) itX.next();
            idxI++; // increment the index

            // get the corresponding network
            MetabolicNetwork startNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(startNetworkSeq));
            if (startNetwork == null) {
                startNetwork = refNetwork.getNetwork(startNetworkSeq);
                lstNetworks.put(strNetwork, startNetwork);
            }
            // get a local copy of directed reactions
            ArrayList lstDirectedReactionsStart = startNetwork.getDirectedReactions();
            Object[] arrDirectedReactionsStart = startNetwork.getDirectedReactions().toArray();

            int idxJ = idxI;
            List partialList = networkList.subList(idxI + 1, nwCount);
            Iterator itY = partialList.iterator();
            while (itY.hasNext()) {
                // get the next network sequence
                Byte[] endNetworkSeq = (Byte[]) itY.next();
                idxJ++; // increment the index

                // get the corresponding network
                MetabolicNetwork endNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(endNetworkSeq));
                if (endNetwork == null) {
                    endNetwork = refNetwork.getNetwork(endNetworkSeq);
                    lstNetworks.put(strNetwork, endNetwork);
                }
                // get a local copy of directed reactions
                ArrayList lstDirectedReactionsEnd = endNetwork.getDirectedReactions();

                ArrayList lstDifferences = MetabolicNetwork.findDifferences(startNetwork, endNetwork);
                if (lstDifferences.size() == 1) {// the networks differ by only one edge
                    int rxnIdxStart = lstDirectedReactionsStart.indexOf(lstDifferences.get(0));
                    int rxnIdxEnd = lstDirectedReactionsEnd.indexOf(lstDifferences.get(0));
                    rateMatrix[idxI][idxJ] = calculateHyperedgeChangeRate(startNetwork, refNetwork, (DirectedReaction) arrDirectedReactionsStart[rxnIdxStart], rxnIdxStart, evolutionModel, insRateQ, delRateQ, dependenceProbability);
//                    rateMatrix[idxI][idxJ] = calculateHyperedgeChangeRate(startNetwork, refNetwork, (DirectedReaction) lstDirectedReactionsStart.get(rxnIdxStart), rxnIdxStart, evolutionModel, insRateQ, delRateQ);
                    rateMatrix[idxJ][idxI] = calculateHyperedgeChangeRate(endNetwork, refNetwork, (DirectedReaction) lstDirectedReactionsEnd.get(rxnIdxEnd), rxnIdxEnd, evolutionModel, insRateQ, delRateQ, dependenceProbability);
                } // networks differ by single edge
            } // end for j
        } // end for i

        // set the diagonal elements
        for (int i = 0; i < rateMatrix.length; i++) {
            double totRate = Utilities.sum(rateMatrix[i]);
            rateMatrix[i][i] = -totRate;
        }

        return rateMatrix;
    }

    public ArrayList enumerateAllNetworks(int edgeCount) {

        int nwCount = (int) Math.pow(2, edgeCount);
        ArrayList lstAllNetworks = new ArrayList(nwCount);
        for (int i = 0; i < nwCount; i++) {
            // get the binary string corresponding to the number i
            String strBinary = Utilities.padString(Integer.toBinaryString(i), edgeCount, '0', true);
            // convert to char array
            char[] strBinaryArr = strBinary.toCharArray();

            // construct the network
            Byte[] theNetwork = new Byte[edgeCount];
            for (int j = 0; j < strBinaryArr.length; j++)
                theNetwork[j] = strBinaryArr[strBinaryArr.length - j - 1] == '1' ? MetabolicNetwork.SEQ_ENTRY_PRESENT : MetabolicNetwork.SEQ_ENTRY_ABSENT;
//                theNetwork[j] = strBinaryArr[j] == '1' ? MetabolicNetwork.SEQ_ENTRY_PRESENT : MetabolicNetwork.SEQ_ENTRY_ABSENT;

            // add to the list
            lstAllNetworks.add(theNetwork);
        //System.out.println(Utilities.toString(theNetwork));
        }

        return lstAllNetworks;
    }

    public ArrayList enumerateAllNetworks(int edgeCount, Byte[] unalterableEdges) {
        // coreEdges: 1s for edges that are core
        // prohibEdges: 1s for edges that are prohib
        // unalterableEdges: 0s for being absent, 1s for being present and -1 for being alterable

        int nwCount = (int) Math.pow(2, edgeCount);
        ArrayList lstAllNetworks = new ArrayList(nwCount);
        for (int i = 0; i < nwCount; i++) {
            // get the binary string corresponding to the number i
            String strBinary = Utilities.padString(Integer.toBinaryString(i), edgeCount, '0', true);
            // convert to char array
            char[] strBinaryArr = strBinary.toCharArray();

            int idx = 0;
            char[] networkSeq = new char[unalterableEdges.length];
            for (int j = 0; j < networkSeq.length; j++) {
                // if its an alterable edge then take its status from binary array
                // otherwise mark it as absent
                if (unalterableEdges[j].equals(new Byte("-1")))
                    networkSeq[j] = strBinaryArr[idx++];
                else
                    networkSeq[j] = '0';
            }

            // construct the network
            Byte[] theNetwork = new Byte[unalterableEdges.length];
            boolean legalNetwork = true;
            for (int j = 0; j < networkSeq.length; j++) {
                theNetwork[j] = networkSeq[networkSeq.length - j - 1] == '1' ? MetabolicNetwork.SEQ_ENTRY_PRESENT : MetabolicNetwork.SEQ_ENTRY_ABSENT;

            }

            // add to the list
            if (legalNetwork)
                lstAllNetworks.add(theNetwork);
        }

        lstAllNetworks.trimToSize();
        return lstAllNetworks;
    }

    public ArrayList enumerateNetworks(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int extraEvents) {

        if (extraEvents < 0)
            return null;

        // find the differences between the two networks
        ArrayList lstDifferences = MetabolicNetwork.findDifferences(startNetwork, endNetwork);
        // number of differences
        int nDifferences = lstDifferences.size();

        // get the edges that can be altered (all edges - core edges - prohib edges)
        ArrayList alterableEdges = startNetwork.getDirectedReactions();
        alterableEdges.removeAll(coreEdges);
        alterableEdges.removeAll(prohibEdges);
        // number of alterable edges
        int nAlterableEdges = alterableEdges.size();

        // Extra Edges (alterable edges - different eges) .. these edges will result in new networks
        ArrayList extraEdges = (ArrayList) alterableEdges.clone();
        extraEdges.removeAll(lstDifferences);

        // get the total number of networks to be seen
        int networkCountShortest = getNetworkCount(nAlterableEdges, nDifferences, 0);
        int networkCount = getNetworkCount(nAlterableEdges, nDifferences, extraEvents);

        // Data structure to hold the networks
        ArrayList networkList = new ArrayList(networkCount);

        // Data structure to hold the networks on shortest path
        ArrayList networkListShortest = new ArrayList(networkCountShortest);
        // add the starting and final networks
        networkListShortest.add(startNetwork.getReactionSequence());
        networkListShortest.add(endNetwork.getReactionSequence());
        // get the networks on shortest path
        networksOnShortestPath(startNetwork, endNetwork, lstDifferences, networkListShortest);


        // add the networks found on the shortest path to the list
        networkList.addAll(networkListShortest);
        // process each of these networks
        Iterator it = networkListShortest.iterator();
        while (it.hasNext()) {
            MetabolicNetwork theNetwork = startNetwork.getNetwork((Byte[]) it.next());

            networksUsingExtraEvents(theNetwork, endNetwork, extraEdges, networkList, extraEvents);
        }
        return networkList;
    }

    private void networksOnShortestPath(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, ArrayList lstEdges, ArrayList networkList) {

        ArrayList alterableEdges = (ArrayList) lstEdges.clone();
        while (alterableEdges.size() > 0) {
            MetabolicNetwork newNetwork = startNetwork.updateNetwork((DirectedReaction) alterableEdges.get(0));

            if (MetabolicNetwork.findDifferences(newNetwork, endNetwork).size() != 0)
                networkList.add(newNetwork.getReactionSequence());

            alterableEdges.remove(0);
            // recursive call for remaining edges
            networksOnShortestPath(newNetwork, endNetwork, alterableEdges, networkList);
        }
    }

    private void networksUsingExtraEvents(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, ArrayList lstEdges, ArrayList networkList, int extraEvents) {

        if (extraEvents == 0)
            return;
        ArrayList alterableEdges = (ArrayList) lstEdges.clone();
        while (alterableEdges.size() > 0) {
            MetabolicNetwork newNetwork = startNetwork.updateNetwork((DirectedReaction) alterableEdges.get(0));

            if (MetabolicNetwork.findDifferences(newNetwork, endNetwork).size() != 0)
                networkList.add(newNetwork.getReactionSequence());

            alterableEdges.remove(0);
            // recursive call for remaining edges
            networksUsingExtraEvents(newNetwork, endNetwork, alterableEdges, networkList, extraEvents - 1);
        }
    }

    public int getNetworkCount(int edgeCount, int nDifferences, int extraEvents) {

        // shortest path
        if (extraEvents == 0)
            return (int) Math.pow(2, nDifferences);

        long numerator = 1;
        int nAlterableEdges = edgeCount - nDifferences;
        for (int i = nAlterableEdges; i > nAlterableEdges - extraEvents; i--)
            numerator *= i;
        long denominator = 1;
        for (int i = 1; i <= extraEvents; i++)
            denominator *= i;
        int factor = (int) (numerator / denominator);

        return factor * (int) Math.pow(2, nDifferences) + getNetworkCount(edgeCount, nDifferences, extraEvents - 1);
    }

//    public double approximateEquilibriumProbability(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork, 
//            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, double insRate, double delRate) {
//
//        // get the change rates for the edges
//        double[] edgeRates = calculateHyperedgeChangeRates(theNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate);
//        
////        // Create a complement network
////        Byte[] rxnSeq = theNetwork.getReactionSequence();
////        Byte[] rxnSeqComplement = new Byte[rxnSeq.length];
////        for (int i = 0; i < rxnSeq.length; i++)
////            rxnSeqComplement[i] = (rxnSeq[i] == MetabolicNetwork.SEQ_ENTRY_PRESENT ? MetabolicNetwork.SEQ_ENTRY_ABSENT : MetabolicNetwork.SEQ_ENTRY_PRESENT);
////        
////        MetabolicNetwork complementNetwork = refNetwork.getNetwork(rxnSeqComplement);
////        // get the change rates for the edges in the complement network
////        double[] edgeRatesComplement = calculateHyperedgeChangeRates(complementNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate);
//
//        // create the unity matrix
//        Matrix one = new Matrix(new double[][]{{1.0,1.0},{1.0,1.0}});
//        
//        double logProb = 0.0;
//        // process all edges
//        int idx = -1;
//        Iterator it = theNetwork.getDirectedReactions().iterator();
//        while (it.hasNext()) {
//            // get the next reaction
//            DirectedReaction reaction = (DirectedReaction)it.next();
//            // increment the index
//            idx++;
//
//            //Byte[] rxnSeqComplement = (Byte[]) rxnSeq.clone();
//            //rxnSeqComplement[idx] = (rxnSeq[idx] == MetabolicNetwork.SEQ_ENTRY_PRESENT ? MetabolicNetwork.SEQ_ENTRY_ABSENT : MetabolicNetwork.SEQ_ENTRY_PRESENT);
//        
//            MetabolicNetwork complementNetwork = theNetwork.updateNetwork(reaction);
//            double edgeRateComplement = calculateHyperedgeChangeRate(complementNetwork, refNetwork, reaction, idx, evolutionModel, insRate, delRate);
//            
//            // set correct rates
//            double insRateQ, delRateQ;
//            if (reaction.isActive()) {
//                //insRateQ = edgeRatesComplement[idx];
//                insRateQ = edgeRateComplement;
//                delRateQ = edgeRates[idx];
//            } else {
//                insRateQ = edgeRates[idx];
//                //delRateQ = edgeRatesComplement[idx];                
//                delRateQ = edgeRateComplement;
//            }
//            // create rate matrix for this edge
//            double[][] rateMatrix = new double[][]{{-insRateQ, insRateQ},{delRateQ, -delRateQ}};
//            
//            Matrix rm = new Matrix(rateMatrix);
//            // solve the system (pi*Q = 0 or Q'*pi = 0)
//            Matrix X = rm.transpose().plus(one).solve(one);
//            // get the equilibrium distribution
//            double[] pi = new double[]{X.get(0, 0), X.get(1, 1)};
//            
//            // add to the log probability
//            logProb += Math.log(reaction.isActive() ? pi[1] : pi[0]);
//        }
//        
//        return logProb;
//    }
    public double approximateEquilibriumProbability(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, double insRate, double delRate, double dependenceProbability) {

        // create the unity matrix
        Matrix one = new Matrix(new double[][]{{1.0, 1.0}, {1.0, 1.0}});

        double logProb = 0.0;
        // process all edges
        int idx = -1;
        Iterator it = theNetwork.getDirectedReactions().iterator();
        while (it.hasNext()) {
            // get the next reaction
            DirectedReaction reaction = (DirectedReaction) it.next();
            // increment the index
            idx++;

            if (coreEdges.contains(reaction) || prohibEdges.contains(reaction))
                continue; // zero probability of changing for core & prohibited edges

            // create rate matrix for this edge
            double[][] edgeRateMatrix = calculateEdgeRateMatrix(theNetwork, refNetwork, reaction, idx, evolutionModel, insRate, delRate, dependenceProbability);

            Matrix rm = new Matrix(edgeRateMatrix);
            // solve the system (pi*Q = 0 or Q'*pi = 0)
            Matrix X = rm.transpose().plus(one).solve(one);
            // get the equilibrium distribution
            double[] pi = new double[]{X.get(0, 0), X.get(1, 1)};

            // add to the log probability
            logProb += Math.log(reaction.isActive() ? pi[1] : pi[0]);
        }

        return logProb;
    }

    public double approximateEquilibriumProbability(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork,
            int evolutionModel, double insRate, double delRate, double dependenceProbability) {

        return approximateEquilibriumProbability(theNetwork, refNetwork, evolutionModel, insRate, delRate, DEFAULT_SUB_NETWORK_SIZE, dependenceProbability);
    }

    public double approximateEquilibriumProbability(MetabolicNetwork theNetwork, MetabolicNetwork refNetwork,
            int evolutionModel, double insRate, double delRate, int subNetworkSize, double dependenceProbability) {

        if (refNetwork.getDirectedReactions() == null)
            refNetwork.setupNetworkSequence();

        // if sub network size is less than zero then use the full network
        if (subNetworkSize <= 0)
            subNetworkSize = refNetwork.getDirectedReactions().size();

        return Math.exp(approximateEquilibriumProbability_NeighbourDependence(theNetwork.getReactionSequence(), refNetwork, evolutionModel, insRate, delRate, subNetworkSize, dependenceProbability));
    }

    public double approximateEquilibriumProbability(Byte[] theNetwork, MetabolicNetwork refNetwork,
            int evolutionModel, double insRate, double delRate, int subNetworkSize, double dependenceProbability) {

        if (refNetwork.getDirectedReactions() == null)
            refNetwork.setupNetworkSequence();

        // if sub network size is less than zero then use the full network
        if (subNetworkSize <= 0)
            subNetworkSize = refNetwork.getDirectedReactions().size();

        return Math.exp(approximateEquilibriumProbability_NeighbourDependence(theNetwork, refNetwork, evolutionModel, insRate, delRate, subNetworkSize, dependenceProbability));
    }

    private double approximateEquilibriumProbability_NeighbourDependence(Byte[] theNetwork, MetabolicNetwork refNetwork,
            int evolutionModel, double insRate, double delRate, int subNetworkSize, double dependenceProbability) {

        int edgeCount = refNetwork.getActiveDirectedReactions().size();
        int totalEdgeCount = refNetwork.getDirectedReactions().size();

        // array to mark hyperedges to be included in the sub networks
        Byte[] refSubNetworkSeq;
        // variables to store new sub networks
        MetabolicNetwork refSubNetwork;
        // array to unalterable edges
        Byte[] unalterableEdges = new Byte[totalEdgeCount];


        boolean subNetworksUsed;
        int nwIdx;
        if (edgeCount == 0)
            return 0.0;
        else if (edgeCount <= subNetworkSize) {

            // no subnetworks used
            subNetworksUsed = false;

            refSubNetworkSeq = refNetwork.getReactionSequence().clone();
            // define edges that are active in the reference network as alterable. 
            // all the rest are unalterable edges
            for (int i = 0; i < totalEdgeCount; i++) {
                if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    unalterableEdges[i] = -1;
                else
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            } // end for all edges
        } else {

            // subnetworks used
            subNetworksUsed = true;

            // initialise
            refSubNetworkSeq = new Byte[totalEdgeCount];
            unalterableEdges = new Byte[totalEdgeCount];
            // initialise
            for (int i = 0; i < totalEdgeCount; i++) {
                refSubNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
                unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            }

            // get the full reference network
            Byte[] allOnes = new Byte[totalEdgeCount];
            for (int i = 0; i < allOnes.length; i++) {
                allOnes[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            }

            // get the neighbours count based on full reference network
            MetabolicNetwork fullRefNetwork = refNetwork.getNetwork(allOnes);
            int[] neighboursCount = fullRefNetwork.getNeighboursCount();
            // get the indices according to the neighbour count
            int[] whichMax = Utilities.whichMaximum(neighboursCount, totalEdgeCount);


//            // remove the ones which are not present in the given reference network
//            for (int i = 0; i < totalEdgeCount; i++) {
//                if (refNetwork.getReactionSequence()[i].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT))
//                    neighboursCount[i] = Integer.MIN_VALUE;
//            } // end for all edges

            Object[] arrDirectedReactionsRef = refNetwork.getDirectedReactions().toArray();
//            ArrayList lstDirectedReactionsRef = refNetwork.getDirectedReactions();
            int currEdgeCount = 0;
            while (currEdgeCount < subNetworkSize) {
                // create an array with rank of each hyperedge based on neighbours count
                int[] ranks = new int[totalEdgeCount];
                for (int i = 0; i < totalEdgeCount; i++)
                    ranks[i] = Integer.MAX_VALUE;
                int rank = 0;
                for (int i = 0; i < whichMax.length; i++) {
                    if (!refNetwork.getReactionSequence()[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT) && !refSubNetworkSeq[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                        ranks[whichMax[i]] = rank++;
                }
                // get the reaction having highest number of neighbours
                int rxnIdx = 0;
                for (int i = 0; i < ranks.length; i++) {
                    if (ranks[i] == 0) {
                        rxnIdx = i;
                        break;
                    }
                }
//                DirectedReaction reaction = (DirectedReaction) lstDirectedReactionsRef.get(rxnIdx);
                DirectedReaction reaction = (DirectedReaction) arrDirectedReactionsRef[rxnIdx];
                // get the neighbours of this reaction in the full reference network
                HashMap neighbours = fullRefNetwork.findReactionNeighbours(reaction, rxnIdx);
                // mark the reaction in sub network
                refSubNetworkSeq[rxnIdx] = 1;
                unalterableEdges[rxnIdx] = -1;
                currEdgeCount++;

                Iterator it = neighbours.keySet().iterator();
                while (currEdgeCount < subNetworkSize && it.hasNext()) {
                    // get the reaction having highest number of neighbours
                    int neighbourIdx = (Integer) it.next();

                    if (refNetwork.getReactionSequence()[neighbourIdx].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT))
                        continue;

                    boolean includeNeighbour = true;
                    if (neighbours.size() > subNetworkSize - 1 && ranks[neighbourIdx] > subNetworkSize - 1)
                        includeNeighbour = false;


                    // mark this edge in the reference and current network if this neighbour is to be included
                    if (includeNeighbour) {
                        refSubNetworkSeq[neighbourIdx] = 1;
                        unalterableEdges[neighbourIdx] = -1;
                        currEdgeCount++;
                    } else {
                        refSubNetworkSeq[neighbourIdx] = 0;
                    }

                } // end while (it.hasNext())
            }
        } // end if edgeCount < MAX_EDGE_COUNT

        // get the reference network based on edges being used for equilibrium probability calculation
        refSubNetwork = refNetwork.getNetwork(refSubNetworkSeq);

        // network index
        StringBuffer strTheNetwork = new StringBuffer();
        for (int i = 0; i < refSubNetworkSeq.length; i++) {
            if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                if (theNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    strTheNetwork.append("1");
                else
                    strTheNetwork.append("0");
            }
        }
        nwIdx = Integer.parseInt(new StringBuffer(strTheNetwork).reverse().toString(), 2);

        // create a rate matrix based on these reactions
        ArrayList lstAllNetworks = enumerateAllNetworks(Math.min(edgeCount, subNetworkSize), unalterableEdges);
        // remove and use recursion
        double[][] rateMatrix = generateRateMatrix(lstAllNetworks, refSubNetwork, evolutionModel, insRate, delRate, dependenceProbability);

        // matrix of 1s
        double[][] unityMatrix = new double[rateMatrix.length][rateMatrix.length];
        for (int i = 0; i < unityMatrix.length; i++) {
            for (int j = 0; j < unityMatrix.length; j++) {
                unityMatrix[i][j] = 1.0;
            }
        }

        Matrix rm = new Matrix(rateMatrix);
        Matrix one = new Matrix(unityMatrix);
        // solve the system (pi*Q = 0 or Q'*pi = 0)
        Matrix X = rm.transpose().plus(one).solve(one);
        double logProb = Math.log(X.get(nwIdx, nwIdx));

        // make the recursive call if we used a sub network
        if (subNetworksUsed) {
            // mark the edges just used as absent in the reference network
            Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
            for (int i = 0; i < newRefNetworkSeq.length; i++) {
                if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            }

            MetabolicNetwork newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
            logProb += approximateEquilibriumProbability_NeighbourDependence(theNetwork, newRefNetwork, evolutionModel, insRate, delRate, subNetworkSize, dependenceProbability);
        }

        return logProb;
    }

// Modified to use hashing to avoid calls to MetabolicNetwork.getNetwork()
//    public double approximateTransitionProbability(MetabolicNetwork startNetwork,
//            MetabolicNetwork endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
//            Byte[] iProhibEdges, int evolutionModel, double insRate, double delRate,
//            double evolTime, int subNetworkSize) {
//
//        if (refNetwork.getDirectedReactions() == null)
//            refNetwork.setupNetworkSequence();
//
//        // if sub network size is less than zero then use the full network
//        if (subNetworkSize <= 0)
//            subNetworkSize = refNetwork.getDirectedReactions().size();
//
//        // mark core and prohibited edges as absent from the reference network ... 
//        // we dont want to use them when generating rate matrices
//        Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
//        for (int i = 0; i < newRefNetworkSeq.length; i++) {
//            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//        }
//        MetabolicNetwork newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
//
//        return Math.exp(approximateTransitionProbabilityByNeighbourhood(startNetwork.getReactionSequence(),
//                endNetwork.getReactionSequence(), newRefNetwork, iCoreEdges, iProhibEdges, 
//                evolutionModel, insRate, delRate, evolTime, subNetworkSize));
//    }
    public double approximateTransitionProbability(MetabolicNetwork startNetwork,
            MetabolicNetwork endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
            Byte[] iProhibEdges, int evolutionModel, double insRate, double delRate,
            double evolTime, int subNetworkSize, double dependenceProbability) {

        String strNetwork;

        if (refNetwork.getDirectedReactions() == null)
            refNetwork.setupNetworkSequence();

        // if sub network size is less than zero then use the full network
        if (subNetworkSize <= 0)
            subNetworkSize = refNetwork.getDirectedReactions().size();

        // mark core and prohibited edges as absent from the reference network ... 
        // we dont want to use them when generating rate matrices
        Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
        for (int i = 0; i < newRefNetworkSeq.length; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        MetabolicNetwork newRefNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(newRefNetworkSeq));
        if (newRefNetwork == null) {
            newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
            lstNetworks.put(strNetwork, newRefNetwork);
        }

        return Math.exp(approximateTransitionProbabilityByNeighbourhood(startNetwork.getReactionSequence(),
                endNetwork.getReactionSequence(), newRefNetwork, iCoreEdges, iProhibEdges,
                evolutionModel, insRate, delRate, evolTime, subNetworkSize, dependenceProbability));
    }

// Modified to use hashing to avoid calls to MetabolicNetwork.getNetwork()
//    public double approximateTransitionProbability(Byte[] startNetwork,
//            Byte[] endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
//            Byte[] iProhibEdges, int evolutionModel, double insRate, double delRate,
//            double evolTime, int subNetworkSize) {
//
//        if (refNetwork.getDirectedReactions() == null)
//            refNetwork.setupNetworkSequence();
//
//        // if sub network size is less than zero then use the full network
//        if (subNetworkSize <= 0)
//            subNetworkSize = refNetwork.getDirectedReactions().size();
//
//        // mark core and prohibited edges as absent from the reference network ... 
//        // we dont want to use them when generating rate matrices
//        Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
//        for (int i = 0; i < newRefNetworkSeq.length; i++) {
//            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//        }
//        MetabolicNetwork newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
//
//        return Math.exp(approximateTransitionProbabilityByNeighbourhood(startNetwork,
//                endNetwork, newRefNetwork, iCoreEdges, iProhibEdges, 
//                evolutionModel, insRate, delRate, evolTime, subNetworkSize));
//    }
    public double approximateTransitionProbability(Byte[] startNetwork,
            Byte[] endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
            Byte[] iProhibEdges, int evolutionModel, double insRate, double delRate,
            double evolTime, int subNetworkSize, double dependenceProbability) {

        String strNetwork;

        if (refNetwork.getDirectedReactions() == null)
            refNetwork.setupNetworkSequence();

        // if sub network size is less than zero then use the full network
        if (subNetworkSize <= 0)
            subNetworkSize = refNetwork.getDirectedReactions().size();

        // mark core and prohibited edges as absent from the reference network ... 
        // we dont want to use them when generating rate matrices
        Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
        for (int i = 0; i < newRefNetworkSeq.length; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        MetabolicNetwork newRefNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(newRefNetworkSeq));
        if (newRefNetwork == null) {
            newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
            lstNetworks.put(strNetwork, newRefNetwork);
        }

        return Math.exp(approximateTransitionProbabilityByNeighbourhood(startNetwork,
                endNetwork, newRefNetwork, iCoreEdges, iProhibEdges,
                evolutionModel, insRate, delRate, evolTime, subNetworkSize, dependenceProbability));
    }

// Modified to use hashing to avoid calls to MetabolicNetwork.getNetwork()
//    private double approximateTransitionProbabilityByNeighbourhood(Byte[] startNetwork,
//            Byte[] endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges, Byte[] iProhibEdges,
//            int evolutionModel, double insRate, double delRate, double evolTime, int subNetworkSize) {
//
//        int edgeCount = refNetwork.getActiveDirectedReactions().size();
//        int totalEdgeCount = refNetwork.getDirectedReactions().size();
//
//        // array to mark hyperedges to be included in the sub networks
//        Byte[] refSubNetworkSeq;
//        // variables to store new sub networks
//        MetabolicNetwork refSubNetwork;
//        // array to unalterable edges
//        Byte[] unalterableEdges = new Byte[totalEdgeCount];
//
//
//        boolean subNetworksUsed;
//        int startIdx, endIdx;
//        if (edgeCount == 0)
//            return 0.0;
//        else if (edgeCount <= subNetworkSize) {
//
//            // no subnetworks used
//            subNetworksUsed = false;
//
//            refSubNetworkSeq = refNetwork.getReactionSequence().clone();
//            // define core and prohibited edges.
//            // also define edges that are active in the reference network as alterable.
//            // all the rest are unalterable edges.
//            for (int i = 0; i < totalEdgeCount; i++) {
//                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//                else if (iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//                else if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    unalterableEdges[i] = -1;
//                else
//                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//            } // end for all edges
//        } else {
//
//            // subnetworks used
//            subNetworksUsed = true;
//
//            // initialise
//            refSubNetworkSeq = new Byte[totalEdgeCount];
//            unalterableEdges = new Byte[totalEdgeCount];
//            // initialise
//            for (int i = 0; i < totalEdgeCount; i++) {
//                refSubNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//
//                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//                else
//                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//            }
//
//            // get the full reference network
//            Byte[] allOnes = new Byte[totalEdgeCount];
//            for (int i = 0; i < allOnes.length; i++) {
//                allOnes[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//            }
//
//            // get the neighbours count based on full reference network
//            MetabolicNetwork fullRefNetwork = refNetwork.getNetwork(allOnes);
//            int[] neighboursCount = fullRefNetwork.getNeighboursCount();
//            // get the indices according to the neighbour count
//            int[] whichMax = Utilities.whichMaximum(neighboursCount, totalEdgeCount);
//
//
//            Object[] arrDirectedReactionsRef = refNetwork.getDirectedReactions().toArray();
////            ArrayList lstDirectedReactionsRef = refNetwork.getDirectedReactions();
//            int currEdgeCount = 0;
//            while (currEdgeCount < subNetworkSize) {
//                // create an array with rank of each hyperedge based on neighbours count
//                int[] ranks = new int[totalEdgeCount];
//                for (int i = 0; i < totalEdgeCount; i++)
//                    ranks[i] = Integer.MAX_VALUE;
//                int rank = 0;
//                for (int i = 0; i < whichMax.length; i++) {
//                    if (!refNetwork.getReactionSequence()[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT) && !refSubNetworkSeq[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                        ranks[whichMax[i]] = rank++;
//                }
//                // get the reaction having highest number of neighbours
//                int rxnIdx = 0;
//                for (int i = 0; i < ranks.length; i++) {
//                    if (ranks[i] == 0) {
//                        rxnIdx = i;
//                        break;
//                    }
//                }
//                DirectedReaction reaction = (DirectedReaction) arrDirectedReactionsRef[rxnIdx];
////                DirectedReaction reaction = (DirectedReaction) lstDirectedReactionsRef.get(rxnIdx);
//                // get the neighbours of this reaction in the full reference network
//                HashMap neighbours = fullRefNetwork.findReactionNeighbours(reaction, rxnIdx);
//                // mark the reaction in sub network
//                refSubNetworkSeq[rxnIdx] = 1;
//                unalterableEdges[rxnIdx] = -1;
//                currEdgeCount++;
//
//                Iterator it = neighbours.keySet().iterator();
//                while (currEdgeCount < subNetworkSize && it.hasNext()) {
//                    // get the reaction having highest number of neighbours
//                    int neighbourIdx = (Integer) it.next();
//
//                    if (refNetwork.getReactionSequence()[neighbourIdx].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT))
//                        continue;
//
//                    boolean includeNeighbour = true;
//                    if (neighbours.size() > subNetworkSize - 1 && ranks[neighbourIdx] > subNetworkSize - 1)
//                        includeNeighbour = false;
//
//
//                    // mark this edge in the reference and current network if this neighbour is to be included
//                    if (includeNeighbour) {
//                        refSubNetworkSeq[neighbourIdx] = 1;
//                        unalterableEdges[neighbourIdx] = -1;
//                        currEdgeCount++;
//                    } else {
//                        refSubNetworkSeq[neighbourIdx] = 0;
//                    }
//
//                } // end while (it.hasNext())
//            }
//        } // end if edgeCount < MAX_EDGE_COUNT
//
//        // network indices
//        StringBuffer strStartNetwork = new StringBuffer();
//        StringBuffer strEndNetwork = new StringBuffer();
//        for (int i = 0; i < refSubNetworkSeq.length; i++) {
//            if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
//                // start network
//                if (startNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    strStartNetwork.append("1");
//                else
//                    strStartNetwork.append("0");
//
//                // end network
//                if (endNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    strEndNetwork.append("1");
//                else
//                    strEndNetwork.append("0");
//            }
//        }
//        startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
//        endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);
//
//        // mark the core and prohibited edges as present 
//        for (int i = 0; i < refSubNetworkSeq.length; i++) {
//            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                refSubNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//        }
//        // get the reference network based on edges being used for this subnetwork 
//        refSubNetwork = refNetwork.getNetwork(refSubNetworkSeq);
//        
//        // create a rate matrix based on these reactions
//        ArrayList lstNetworks = enumerateAllNetworks(Math.min(edgeCount, subNetworkSize), unalterableEdges);
//        double[][] rateMatrix = generateRateMatrix(lstNetworks, refSubNetwork, evolutionModel, insRate, delRate);
////        BigDecimal[][] rateMatrix = generateRateMatrix(lstNetworks, refSubNetwork, evolutionModel, new BigDecimal(insRate), new BigDecimal(delRate));
//
//        // exponentiate to get the transition probabilities
//        double[][] transProbs = Utilities.exponentiate(rateMatrix, evolTime, 12);
////        BigDecimal[][] transProbs = Utilities.exponentiate(rateMatrix, evolTime, 12, 128);
//        // get the log probability from start to end network
//        double logProb = Math.log(transProbs[startIdx][endIdx]);
////        double logProb = Math.log(transProbs[startIdx][endIdx].doubleValue());
//
////        // testing purpose
////        System.out.println(Utilities.toString(rateMatrix));
////        System.out.println();
////        System.out.println(Utilities.toString(transProbs));
//
//        // make the recursive call if we used a sub network
//        if (subNetworksUsed) {
//            // mark the edges just used as absent in the reference network
//            // mark the core and prohibited edges as absent
//            Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
//            for (int i = 0; i < newRefNetworkSeq.length; i++) {
//                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//                else if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
//                    newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//            }
//
//            MetabolicNetwork newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
//            logProb += approximateTransitionProbabilityByNeighbourhood(startNetwork, endNetwork,
//                    newRefNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
//                    evolTime, subNetworkSize);
//        }
//
//        return logProb;
//    }
    private double approximateTransitionProbabilityByNeighbourhood(Byte[] startNetwork,
            Byte[] endNetwork, MetabolicNetwork refNetwork, Byte[] iCoreEdges, Byte[] iProhibEdges,
            int evolutionModel, double insRate, double delRate, double evolTime, int subNetworkSize,
            double dependenceProbability) {

        String strNetwork;
        int edgeCount = refNetwork.getActiveDirectedReactions().size();
        int totalEdgeCount = refNetwork.getDirectedReactions().size();

        // array to mark hyperedges to be included in the sub networks
        Byte[] refSubNetworkSeq;
        // variables to store new sub networks
        MetabolicNetwork refSubNetwork;
        // array to store unalterable edges
        Byte[] unalterableEdges = new Byte[totalEdgeCount];


        boolean subNetworksUsed;
        int startIdx, endIdx;
        if (edgeCount == 0)
            return 0.0;
        else if (edgeCount <= subNetworkSize) {

            // no subnetworks used
            subNetworksUsed = false;

            refSubNetworkSeq = refNetwork.getReactionSequence().clone();
            // define core and prohibited edges.
            // also define edges that are active in the reference network as alterable.
            // all the rest are unalterable edges.
            for (int i = 0; i < totalEdgeCount; i++) {
                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
                else if (iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
                else if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    unalterableEdges[i] = -1;
                else
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            } // end for all edges
        } else {

            // subnetworks used
            subNetworksUsed = true;

            // initialise
            refSubNetworkSeq = new Byte[totalEdgeCount];
            unalterableEdges = new Byte[totalEdgeCount];
            // initialise
            for (int i = 0; i < totalEdgeCount; i++) {
                refSubNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;

                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
                else
                    unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            }

            // get the full reference network
            Byte[] allOnes = new Byte[totalEdgeCount];
            for (int i = 0; i < allOnes.length; i++) {
                allOnes[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            }

            // get the neighbours count based on full reference network
            MetabolicNetwork fullRefNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(allOnes));
            if (fullRefNetwork == null) {
                fullRefNetwork = refNetwork.getNetwork(allOnes);
                lstNetworks.put(strNetwork, fullRefNetwork);
            }

            int[] neighboursCount = fullRefNetwork.getNeighboursCount();
            // get the indices according to the neighbour count
            int[] whichMax = Utilities.whichMaximum(neighboursCount, totalEdgeCount);


            Object[] arrDirectedReactionsRef = refNetwork.getDirectedReactions().toArray();
//            ArrayList lstDirectedReactionsRef = refNetwork.getDirectedReactions();
            int currEdgeCount = 0;
            while (currEdgeCount < subNetworkSize) {
                // create an array with rank of each hyperedge based on neighbours count
                int[] ranks = new int[totalEdgeCount];
                for (int i = 0; i < totalEdgeCount; i++)
                    ranks[i] = Integer.MAX_VALUE;
                int rank = 0;
                for (int i = 0; i < whichMax.length; i++) {
                    if (!refNetwork.getReactionSequence()[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT) && !refSubNetworkSeq[whichMax[i]].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                        ranks[whichMax[i]] = rank++;
                }
                // get the reaction having highest number of neighbours
                int rxnIdx = 0;
                for (int i = 0; i < ranks.length; i++) {
                    if (ranks[i] == 0) {
                        rxnIdx = i;
                        break;
                    }
                }
                DirectedReaction reaction = (DirectedReaction) arrDirectedReactionsRef[rxnIdx];
//                DirectedReaction reaction = (DirectedReaction) lstDirectedReactionsRef.get(rxnIdx);
                // get the neighbours of this reaction in the full reference network
                HashMap neighbours = fullRefNetwork.findReactionNeighbours(reaction, rxnIdx);
                // mark the reaction in sub network
                refSubNetworkSeq[rxnIdx] = 1;
                unalterableEdges[rxnIdx] = -1;
                currEdgeCount++;

                Iterator it = neighbours.keySet().iterator();
                while (currEdgeCount < subNetworkSize && it.hasNext()) {
                    // get the reaction having highest number of neighbours
                    int neighbourIdx = (Integer) it.next();

                    if (refNetwork.getReactionSequence()[neighbourIdx].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT))
                        continue;

                    boolean includeNeighbour = true;
                    if (neighbours.size() > subNetworkSize - 1 && ranks[neighbourIdx] > subNetworkSize - 1)
                        includeNeighbour = false;


                    // mark this edge in the reference and current network if this neighbour is to be included
                    if (includeNeighbour) {
                        refSubNetworkSeq[neighbourIdx] = 1;
                        unalterableEdges[neighbourIdx] = -1;
                        currEdgeCount++;
                    } else {
                        refSubNetworkSeq[neighbourIdx] = 0;
                    }

                } // end while (it.hasNext())
            }
        } // end if edgeCount < MAX_EDGE_COUNT

        // network indices
        StringBuffer strStartNetwork = new StringBuffer();
        StringBuffer strEndNetwork = new StringBuffer();
        for (int i = 0; i < refSubNetworkSeq.length; i++) {
            if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                // start network
                if (startNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    strStartNetwork.append("1");
                else
                    strStartNetwork.append("0");

                // end network
                if (endNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    strEndNetwork.append("1");
                else
                    strEndNetwork.append("0");
            }
        }
        startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
        endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);

        // mark the core and prohibited edges as present 
        for (int i = 0; i < refSubNetworkSeq.length; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                refSubNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        // get the reference network based on edges being used for this subnetwork 
        refSubNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(refSubNetworkSeq));
        if (refSubNetwork == null) {
            refSubNetwork = refNetwork.getNetwork(refSubNetworkSeq);
            lstNetworks.put(strNetwork, refSubNetwork);
        }

        // create a rate matrix based on these reactions
        ArrayList lstAllNetworks = enumerateAllNetworks(Math.min(edgeCount, subNetworkSize), unalterableEdges);
        double[][] rateMatrix = generateRateMatrix(lstAllNetworks, refSubNetwork, evolutionModel, insRate, delRate, dependenceProbability);
//        BigDecimal[][] rateMatrix = generateRateMatrix(lstNetworks, refSubNetwork, evolutionModel, new BigDecimal(insRate), new BigDecimal(delRate));

        // exponentiate to get the transition probabilities
        double[][] transProbs = Utilities.exponentiate(rateMatrix, evolTime, 12);
//        BigDecimal[][] transProbs = Utilities.exponentiate(rateMatrix, evolTime, 12, 128);
        // get the log probability from start to end network
        double logProb = Math.log(transProbs[startIdx][endIdx]);
//        double logProb = Math.log(transProbs[startIdx][endIdx].doubleValue());

//        // testing purpose
//        System.out.println(Utilities.toString(rateMatrix));
//        System.out.println();
//        System.out.println(Utilities.toString(transProbs));

        // make the recursive call if we used a sub network
        if (subNetworksUsed) {
            // mark the edges just used as absent in the reference network
            // mark the core and prohibited edges as absent
            Byte[] newRefNetworkSeq = (Byte[]) refNetwork.getReactionSequence().clone();
            for (int i = 0; i < newRefNetworkSeq.length; i++) {
                if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
                else if (refSubNetworkSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                    newRefNetworkSeq[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            }

            MetabolicNetwork newRefNetwork = (MetabolicNetwork) lstNetworks.get(strNetwork = getNetworkString(newRefNetworkSeq));
            if (newRefNetwork == null) {
                newRefNetwork = refNetwork.getNetwork(newRefNetworkSeq);
                lstNetworks.put(strNetwork, newRefNetwork);
            }
            logProb += approximateTransitionProbabilityByNeighbourhood(startNetwork, endNetwork,
                    newRefNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
                    evolTime, subNetworkSize, dependenceProbability);
        }

        return logProb;
    }

    public double[][] calculateTransitionProbabilityMatrix(double[][] rateMatrix, double evolTime) {

        // exponentiate rate matrix to get the transition probabilities
        return Utilities.exponentiate(rateMatrix, evolTime, 10);
    }

    public double[] calculateEquilibriumProbabilityVector(double[][] rateMatrix) {

        // matrix of 1s
        double[][] unityMatrix = new double[rateMatrix.length][rateMatrix.length];
        for (int i = 0; i < unityMatrix.length; i++) {
            for (int j = 0; j < unityMatrix.length; j++) {
                unityMatrix[i][j] = 1.0;
            }
        }

        Matrix rm = new Matrix(rateMatrix);
        Matrix one = new Matrix(unityMatrix);
        // solve the system (pi*Q = 0 or Q'*pi = 0)
        Matrix X = rm.transpose().plus(one).solve(one);

        // get the equilibrium distribution
        double[] pi = new double[rateMatrix.length];
        for (int i = 0; i < pi.length; i++) {
            pi[i] = X.get(i, i);
        }

        return pi;
    }

    public double[] calculateEquilibriumProbability(double[][] rateMatrix) {

        Matrix rm = new Matrix(rateMatrix);
        // create the unity matrix
        Matrix one = new Matrix(new double[][]{{1.0, 1.0}, {1.0, 1.0}});
        // solve the system (pi*Q = 0 or Q'*pi = 0)
        Matrix X = rm.transpose().plus(one).solve(one);
        // get the equilibrium distribution
        double[] pi = new double[]{X.get(0, 0), X.get(1, 1)};

        return pi;
    }

    private String getNetworkString(Byte[] theNetwork) {
        return Utilities.toString(theNetwork, "");
    }

    static public int getNetworkIndex(Byte[] theNetwork) {
        // network index
        StringBuffer strTheNetwork = new StringBuffer();
        for (int i = 0; i < theNetwork.length; i++) {
            if (theNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                strTheNetwork.append("1");
            else
                strTheNetwork.append("0");
        }
        return Integer.parseInt(new StringBuffer(strTheNetwork).reverse().toString(), 2);

    }

    public double calculateFullLikelihood(PhyloNode phylonode, MetabolicNetwork refNetwork,
            ArrayList lstNetworks, Byte[] parentNetwork, double insRate, double delRate,
            double[][] transProbs, double[] eqProbs) {
        /** this function calculates the likelihood by using full network. 
         *  Working: For given rates,  it first calculates the exponential of the matrix and 
         *           saves it for future use. It then recursively calculates the likelihood of 
         *           the tree by summing over all possible networks
         *  Caution: Not to be used on large networks
         */
        Iterator it;
        if (phylonode.isLeaf()) {
            ArrayList dummyList = new ArrayList();
            dummyList.add(phylonode.getMetabolicNetwork().getReactionSequence());
            it = dummyList.iterator();
        } else
            it = lstNetworks.iterator();

        double treeProb = 0.0;
        while (it.hasNext()) {
            // get the current network
            Byte[] currentNetwork = (Byte[]) it.next();

            Double parentTP;
            if (parentNetwork == null) { // this is the root .. calculate eq. probability
                parentTP = eqProbs[getNetworkIndex(currentNetwork)];
            } else {
                // transition probability from parent to current network
                parentTP = transProbs[getNetworkIndex(parentNetwork)][getNetworkIndex(currentNetwork)];

            } // end if parentNetwork == null

            Double leftTP = 1.0, rightTP = 1.0;
            if (phylonode.getLeftSon() != null)
                // calculate transition probability from current network to Left son
                leftTP = calculateFullLikelihood(phylonode.getLeftSon(), refNetwork, lstNetworks,
                        currentNetwork, insRate, delRate, transProbs, eqProbs);
            if (phylonode.getRightSon() != null)
                // calculate transition probability from current network to Right son
                rightTP = calculateFullLikelihood(phylonode.getRightSon(), refNetwork, lstNetworks,
                        currentNetwork, insRate, delRate, transProbs, eqProbs);

            // add the log-values
            treeProb += Math.exp(Math.log(parentTP) + Math.log(leftTP) + Math.log(rightTP));
        } // end while

        return treeProb;
    }
}
