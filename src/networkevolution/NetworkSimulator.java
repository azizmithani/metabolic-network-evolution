/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package networkevolution;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import network.DirectedReaction;
import network.MetabolicNetwork;
import network.Metabolite;
import network.NetworkObject;
import systemobject.Organism;
import rahnumadatabase.Database;
import utilities.Utilities;

/**
 *
 * @author mithani
 */
public class NetworkSimulator {

    static public String DB_DRIVER = "com.mysql.jdbc.Driver";
//    // test environment
//    static public String DB_URL = "jdbc:mysql://localhost:3306/";
//    static public String DB_NAME = "rahnumatest";
//    static public String DB_USER = "rahnumatestdba";
//    static public String DB_PASSWORD = "oveaQue5";

    // production environment
    static public String DB_URL = "jdbc:mysql://peng.stats.ox.ac.uk:3306/";
    static public String DB_NAME = "rahnuma";
    static public String DB_USER = "rahnumadba";
    static public String DB_PASSWORD = "eetha5Ku";

    static public boolean verbose = true;
    /**
     * @param args the command line arguments
     */
//    public static void main(String[] args) {
//
//        String[] orgIds, pathwayCodes;
//
//        orgIds = new String[]{"pae", "pau", "ppu", "ppf", "pst", "psb", "psp", "pfl", "pfo", "pen", "pmy"};
//        //createDistanceMatrix(orgIds);
//
//        orgIds = new String[]{"pae", "pau", "ppu", "ppf", "pst", "psb", "psp", "pfl", "pfo", "pen", "pmy", "xfa", "xft", "xfm", "xcc", "xcb", "xcv"};
//        pathwayCodes = new String[]{"00271"};
//        //getNeighboursCount(orgIds, pathwayCodes);
//        //getDirectedReactionStatus(orgIds, pathwayCodes);
//        //simulateEvolution(orgIds, pathwayCodes);
//        compareEvolutionForDifferentRates(NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, null, pathwayCodes);
//        compareEvolutionForDifferentRates(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, null, pathwayCodes);
//
////        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
////        database.readApplicationParameters();
////
////        // read the kegg pathways from the db
////        String[][] keggPathways = new String[][]{{"00271", "00271"}};
////        String[] pathwayIds = new String[]{keggPathways[0][0]};
////
////        String[][] orgReactions = null;
////        String[][] rxnPathways = database.getReactionPathways(pathwayIds);
////        
////        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
////                "Reference", new ArrayList(), database.getMetabolites(""),
////                database.getReactions(pathwayIds), rxnPathways, orgReactions, keggPathways,
////                null, null, null, null);
////        
////        rxnPathways = null;
////        // first network
////        Organism org = new Organism("pfo", "pfo");
////        ArrayList orgList = new ArrayList();
////        orgList.add(org);
////        orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
////        MetabolicNetwork orgNetwork1 = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
////                org.getName(), orgList, database.getMetabolites(""),
////                database.getReactions(pathwayIds), rxnPathways, orgReactions, keggPathways,
////                null, null, null, null);
////
////        // second network
////        org = new Organism("pst", "pst");
////        orgList = new ArrayList();
////        orgList.add(org);
////        orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
////        MetabolicNetwork orgNetwork2 = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
////                org.getName(), orgList, database.getMetabolites(""),
////                database.getReactions(pathwayIds), rxnPathways, orgReactions, keggPathways,
////                null, null, null, null);
////
////        
////        orgNetwork1.addInactiveReactions(refNetwork.getReactions());
////        orgNetwork2.addInactiveReactions(refNetwork.getReactions());
////        
////        orgNetwork1.removeInactiveMetabolites();
////        orgNetwork2.removeInactiveMetabolites();
////        refNetwork.removeInactiveMetabolites();
////        
////        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
////        String ignoreFile = "ignoreCompoundsTCA.txt";
////        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));       
////        orgNetwork1.getMetabolites().removeAll(ignoreMetabolites);
////        orgNetwork2.getMetabolites().removeAll(ignoreMetabolites);
////        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
////        
////        orgNetwork1.setupNetworkMatrices();
////        orgNetwork2.setupNetworkMatrices();
////        refNetwork.setupNetworkMatrices();
////
////        orgNetwork1.setupNetworkSequence();
////        orgNetwork2.setupNetworkSequence();
////        
////        int nIter = 0;
////        double insRate = 0.5, delRate = 0.3;
////        ArrayList independentEdge = new Controller().simulateNetworkEvolution(orgNetwork1, EVOL_MODEL_INDEPENDENT_EDGE, nIter, insRate, delRate);
////        new Controller().summariseEvolution(orgNetwork1, independentEdge);
////        
////        ArrayList neighbourDependent = new Controller().simulateNetworkEvolution(orgNetwork1, EVOL_MODEL_NEIGHBOUR_DEPENDENT, nIter, insRate, delRate);       
////        new Controller().summariseEvolution(orgNetwork1, neighbourDependent);
////        int[] getNeighboursCount = new Controller().getNeighboursCount(refNetwork);
//
////        Byte[] seq1 = orgNetwork1.getReactionSequence();
////        Byte[] seq2 = orgNetwork2.getReactionSequence();
////
////        System.out.println(Utilities.toString(NetworkObject.getIDs(orgNetwork1.getDirectedReactions())));
////        System.out.println(Utilities.toString(seq1));
////        System.out.println(Utilities.toString(seq2));
//
////                orgNetwork.removeInactiveMetabolites();
////                orgNetwork.setupNetworkMatrices();
////                Object[][] nei = orgNetwork.getNodeEdgeIncidenceMatrix().full(Byte.class);
////                for (int k = 0; k < nei.length; k++) {
////                    Object[] row = nei[k];
////                    for (int l = 0; l < row.length; l++) {
////                        System.out.print((Byte) row[l]);
////                        System.out.print("\t");
////                    }
////                    System.out.println();
////                }
////                System.out.println();
////                Object[][] nb = orgNetwork.getNeighbourhoodMatrix().full(Byte.class);
////                for (int k = 0; k < nb.length; k++) {
////                    Object[] row = nb[k];
////                    for (int l = 0; l < row.length; l++) {
////                        System.out.print((Byte) row[l]);
////                        System.out.print("\t");
////                    }
////                    System.out.println();
////                }
//
//    }

    public ArrayList simulateNetworkEvolution(MetabolicNetwork theNetwork, int model, int nIter, 
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, int updateInterval, double dependenceProbability) {
        
        return simulateNetworkEvolution(theNetwork, model, nIter, 0, coreEdges, prohibEdges, insRate, delRate, updateInterval, dependenceProbability);
    }
    
    public ArrayList simulateNetworkEvolution(MetabolicNetwork theNetwork, int model, int nIter, int nBurning,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, int updateInterval,
            double dependenceProbability) {

        if (nBurning > 0) {
            System.out.println("Burn-in Period Started");
            for (int i = 0; i < nBurning; i++) {
                DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
                theNetwork = theNetwork.updateNetwork(hyperedge);
            }
            System.out.println("Burn-in Period Finished");
        }
        
        ArrayList lstNetworks = new ArrayList(nIter - nBurning + 1);
        lstNetworks.add(theNetwork.getReactionSequence());
        System.out.println("Simulation Started");
        int size = 0;
        for (int i = 0; i < nIter; i++) {
//            System.out.println(i);
            if (updateInterval > 0 && i % updateInterval ==0)
                System.out.println("Iteration: " + Integer.toString(i) + "\t" + Utilities.toString(theNetwork.getReactionSequence(), " "));
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);
            lstNetworks.add(theNetwork.getReactionSequence());
//            System.out.print(Integer.toString(i) + ":\t");
//            System.out.println(Utilities.toString(theNetwork.getReactionSequence(), " "));
            size += theNetwork.getActiveDirectedReactions().size();
        }
        System.out.println("Simulation Finished");
        System.out.println("Average Network Size: " + Double.toString((double)size/(double)(lstNetworks.size()-1)));
        return lstNetworks;
    }

    public ArrayList simulateNetworkEvolutionNeighboursCount(MetabolicNetwork theNetwork, 
            int model, int nIter, int nBurning, ArrayList coreEdges, ArrayList prohibEdges,
            double insRate, double delRate, int updateInterval, double dependenceProbability) {
            // prints number of neighbours for each hyperedge at each iteration instead of networks

        if (nBurning > 0) {
            System.out.println("Burn-in Period Started");
            for (int i = 0; i < nBurning; i++) {
                DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
                theNetwork = theNetwork.updateNetwork(hyperedge);
            }
            System.out.println("Burn-in Period Finished");
        }

        ArrayList lstNetworks = new ArrayList(nIter - nBurning + 1);
        lstNetworks.add(theNetwork.getReactionSequence());
        System.out.println("Simulation Started");
        int size = 0;
        for (int i = 0; i < nIter; i++) {
//            System.out.println(i);
            if (updateInterval > 0 && i % updateInterval ==0)
                System.out.println("Iteration: " + Integer.toString(i) + "\t" + Utilities.toString(theNetwork.getNeighboursCount(), " "));
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);
            lstNetworks.add(theNetwork.getReactionSequence());
//            System.out.print(Integer.toString(i) + ":\t");
//            System.out.println(Utilities.toString(theNetwork.getReactionSequence(), " "));
            size += theNetwork.getActiveDirectedReactions().size();
        }
        System.out.println("Simulation Finished");
        return lstNetworks;
    }

    public ArrayList simulateNetworkEvolutionWithTime(MetabolicNetwork theNetwork, int model, double evolTime,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, double dependenceProbability) {

        ArrayList lstNetworks = new ArrayList();
        ArrayList lstTimes = new ArrayList();
        double[] edgeRates = new double[theNetwork.getReactionSequence().length];
        
        int idx = 0;
        double time = 0.0;
        lstNetworks.add(theNetwork.getReactionSequence());
        lstTimes.add(time);
        while (time < evolTime) {
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, edgeRates, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);

            // get the total rate
            double totalRate = Utilities.sum(edgeRates);
            
            // Simulate the next event time
            double u_time = new Random().nextDouble();
            double dt = - Math.log(u_time)/totalRate; 
    
            // Evolution time so far
            time += dt;

            // add this network only if time is less than
            if (time < evolTime) {
                lstNetworks.add(theNetwork.getReactionSequence());
                lstTimes.add(time);
            }
            if (verbose) {
                System.out.print(Integer.toString(idx++) + ":\t" + Double.toString(time) + "\t");
                System.out.println(Utilities.toString(theNetwork.getReactionSequence()));
            }
            
        }
        
        ArrayList results = new ArrayList();
        results.add(lstNetworks);
        results.add(lstTimes);
        return results;
    }

    public ArrayList simulateNetworkEvolutionWithTime(MetabolicNetwork theNetwork, int model, double evolTime, 
            int nBurning, ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, 
            double dependenceProbability) {

        ArrayList lstNetworks = new ArrayList();
        ArrayList lstTimes = new ArrayList();
        double[] edgeRates = new double[theNetwork.getReactionSequence().length];
        
        System.out.println("Burn-in Period Started");
        for (int i = 0; i < nBurning; i++) {
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);
        }
        System.out.println("Burn-in Period Finished");
        
        int insCount = 0;
        int delCount = 0;
        int idx = 0;
        double time = 0.0;
        lstNetworks.add(theNetwork.getReactionSequence());
        lstTimes.add(time);
        while (time < evolTime) {
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, edgeRates, dependenceProbability);
            if (hyperedge.isActive()) 
                delCount++;
            else
                insCount++;

            //update the network
            theNetwork = theNetwork.updateNetwork(hyperedge);
            
            // get the total rate
            double totalRate = Utilities.sum(edgeRates);
            
            // Simulate the next event time
            double u_time = new Random().nextDouble();
            double dt = - Math.log(u_time)/totalRate; 
    
            // Evolution time so far
            time += dt;

            // add this network only if time is less than
            if (time < evolTime) {
                lstNetworks.add(theNetwork.getReactionSequence());
                lstTimes.add(time);
            }
            if (verbose) {
                System.out.print(Integer.toString(idx++) + ":\t" + Double.toString(time) + "\t");
                System.out.println(Utilities.toString(theNetwork.getReactionSequence()));
            }
            
        }
        
        ArrayList results = new ArrayList();
        results.add(lstNetworks);
        results.add(lstTimes);
        results.add(insCount);
        results.add(delCount);
        return results;
    }

    public ArrayList simulateNetworkEvolutionWithTime(MetabolicNetwork theNetwork, int model, int nIter, int nBurning,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, double dependenceProbability) {

        ArrayList lstNetworks = new ArrayList();
        ArrayList lstTimes = new ArrayList();
        double[] edgeRates = new double[theNetwork.getReactionSequence().length];
        
        System.out.println("Burn-in Period Started");
        for (int i = 0; i < nBurning; i++) {
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);
        }
        System.out.println("Burn-in Period Finished");
        
        int idx = 0;
        double time = 0.0;
        lstNetworks.add(theNetwork.getReactionSequence());
        lstTimes.add(time);
        for (int i = 0; i < nIter; i++) {
            DirectedReaction hyperedge = selectHyperedge(theNetwork, model, coreEdges, prohibEdges, insRate, delRate, edgeRates, dependenceProbability);
            theNetwork = theNetwork.updateNetwork(hyperedge);

            // get the total rate
            double totalRate = Utilities.sum(edgeRates);
            
            // Simulate the next event time
            double u_time = new Random().nextDouble();
            double dt = - Math.log(u_time)/totalRate; 
    
            // Evolution time so far
            time += dt;

            // add to the list
            lstNetworks.add(theNetwork.getReactionSequence());
            lstTimes.add(time);

            if (verbose) {
                System.out.print(Integer.toString(idx++) + ":\t" + Double.toString(time) + "\t");
                System.out.println(Utilities.toString(theNetwork.getReactionSequence()));
            }
            
        }
        
        ArrayList results = new ArrayList();
        results.add(lstNetworks);
        results.add(lstTimes);
        return results;
    }
    
    static public DirectedReaction selectHyperedge(MetabolicNetwork theNetwork, int model,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, 
            double dependenceProbability) {

        NetworkEvolver networkEvolver = new NetworkEvolver();
        double[] edgeProb = networkEvolver.calculateHyperedgeProbabilities(theNetwork, coreEdges, prohibEdges, model, insRate, delRate, dependenceProbability);

        return (DirectedReaction) theNetwork.getDirectedReactions().get(Utilities.randsample(edgeProb));
    }

    static public DirectedReaction selectHyperedge(MetabolicNetwork theNetwork, int evolutionModel,
        ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate, double[] edgeRates,
        double dependenceProbability) {
        
        NetworkEvolver networkEvolver = new NetworkEvolver();
        double[] rates = networkEvolver.calculateHyperedgeChangeRates(theNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);

        // copy the rates
        System.arraycopy(rates, 0, edgeRates, 0, rates.length);
        
        // get the sum of rates
        double totalRate = Utilities.sum(edgeRates);
        // normalise the rates to get probabilities
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

        // select a hyperedge based on edge probabilities
        DirectedReaction hyperedge = (DirectedReaction) theNetwork.getDirectedReactions().get(Utilities.randsample(edgeProbabilities));

        return hyperedge;
    }
    
    static public int[][] getNeighboursCount(String[] orgIds, String[] pathwayIds) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
        String ignoreFile = "ignoreCompoundsTCA.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));


        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // org networks + ref network
        int[][] neighboursCount = new int[refNetwork.getDirectedReactions().size()][orgIds.length + 1];

        // add the neighbours count for reference network
        int[] refNeighboursCount = refNetwork.getNeighboursCount();
        for (int j = 0; j < refNeighboursCount.length; j++) {
            neighboursCount[j][0] = refNeighboursCount[j];
        }

        // add the neighbours count for org networks
        for (int i = 0; i < orgIds.length; i++) {
            String[][] arrOrg = database.getOrganisms(new String[]{orgIds[i]});
            Organism org = new Organism(arrOrg[0][0], arrOrg[0][1]);
            ArrayList orgList = new ArrayList();
            orgList.add(org);

            String[][] orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(new String[]{org.getID()}, pathwayIds);

            MetabolicNetwork orgNetwork = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
                    org.getName(), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            orgNetwork.addInactiveReactions(refNetwork.getReactions());
            orgNetwork.removeInactiveMetabolites();
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);

            orgNetwork.setupNetworkMatrices();
            orgNetwork.setupNetworkSequence();

            int[] orgNeighboursCount = orgNetwork.getNeighboursCount(refNetwork);
            for (int j = 0; j < orgNeighboursCount.length; j++) {
                neighboursCount[j][i + 1] = orgNeighboursCount[j];
            }
        }


        System.out.println("Neighbours Count:");
        System.out.print("ID" + "\t" + "REF" + "\t");
        for (int i = 0; i < orgIds.length; i++)
            System.out.print(orgIds[i] + "\t");
        System.out.println();
        for (int i = 0; i < neighboursCount.length; i++) {
            int[] is = neighboursCount[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < is.length; j++) {
                System.out.print(is[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        return neighboursCount;
    }

    static public int[][] getDirectedReactionStatus(String[] orgIds, String[] pathwayIds) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
        String ignoreFile = "ignoreCompoundsTCA.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        refNetwork.setupNetworkSequence();

        // data of all org networks 
        int[][] reactionStatus = new int[refNetwork.getDirectedReactions().size()][orgIds.length];

        // add the reaction status for org networks
        for (int i = 0; i < orgIds.length; i++) {
            String[][] arrOrg = database.getOrganisms(new String[]{orgIds[i]});
            Organism org = new Organism(arrOrg[0][0], arrOrg[0][1]);
            ArrayList orgList = new ArrayList();
            orgList.add(org);

            String[][] orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(new String[]{org.getID()}, pathwayIds);

            MetabolicNetwork orgNetwork = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
                    org.getName(), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            orgNetwork.addInactiveReactions(refNetwork.getReactions());
            orgNetwork.removeInactiveMetabolites();
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);

            orgNetwork.setupNetworkSequence();

            Byte[] reactionSequence = orgNetwork.getReactionSequence();
            for (int j = 0; j < reactionSequence.length; j++) {
                reactionStatus[j][i] = (reactionSequence[j].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) ? 1 : 0);
            }
        }


        System.out.println("Reaction Status:");
        System.out.print("ID" + "\t");
        for (int i = 0; i < orgIds.length; i++)
            System.out.print(orgIds[i] + "\t");
        System.out.println();
        for (int i = 0; i < reactionStatus.length; i++) {
            int[] is = reactionStatus[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < is.length; j++) {
                System.out.print(is[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        return reactionStatus;
    }

    static public void simulateEvolution(String[] orgIds, String[] pathwayIds, String[] orgIdsCore, String[] orgIdsProhib, int updateInterval) {

        int nIter = 10000;
        double insRate = 0.5, delRate = 0.3;

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        if (orgIds == null)
            orgIds = new String[]{};

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
        String ignoreFile = "ignoreCompoundsTCA.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();


        // get the list of core reactions
        ArrayList coreReactions;
        if (orgIdsCore != null) {
            String[][] arrOrgCore = database.getOrganisms(orgIdsCore);
            ArrayList orgList = Organism.buildOrganismList(arrOrgCore);

            String[][] orgReactions = database.getOrganismReactions(orgIdsCore, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIdsCore, pathwayIds);

            MetabolicNetwork coreNetwork = MetabolicNetwork.buildMetabolicNetwork("core",
                    "Core", orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            // core reactions = reactions present in all orgIdsCore
            coreNetwork.markCoreReactions();
            coreReactions = coreNetwork.getCoreReactions();
        } else
            coreReactions = new ArrayList();

        // get the list of prohibted reactions
        ArrayList prohibReactions;
        if (orgIdsProhib != null) {
            String[][] arrOrgProhib = database.getOrganisms(orgIdsProhib);
            ArrayList orgList = Organism.buildOrganismList(arrOrgProhib);

            String[][] orgReactions = database.getOrganismReactions(orgIdsProhib, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIdsProhib, pathwayIds);

            MetabolicNetwork prohibNetwork = MetabolicNetwork.buildMetabolicNetwork("prohib",
                    "Prohibited", orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            // prohibited reactions = ref n/w - reactions present in atleast one of the org. in orgIdsProhib
            prohibReactions = (ArrayList) refNetwork.getReactions().clone();
            prohibReactions.removeAll(prohibNetwork.getReactions());
        } else
            prohibReactions = new ArrayList();

        int[][] insFreq = new int[refNetwork.getDirectedReactions().size()][orgIds.length + 1]; // number of insertions for each edge
        int[][] delFreq = new int[refNetwork.getDirectedReactions().size()][orgIds.length + 1]; // number of deletions for each edge

        NetworkSimulator networkSimulator = new NetworkSimulator();
        ArrayList lstNetworks = networkSimulator.simulateNetworkEvolution(refNetwork, NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, nIter, new ArrayList(), new ArrayList(), insRate, delRate, updateInterval, 1.0);
        // variables to hold summary data 
        int[] insFreqRef = new int[refNetwork.getDirectedReactions().size()]; // number of insertions for each edge
        int[] delFreqRef = new int[refNetwork.getDirectedReactions().size()]; // number of deletions for each edge
        NetworkEvolver.summariseEvolution(lstNetworks, insFreqRef, delFreqRef, null, null);

        for (int j = 0; j < insFreqRef.length; j++) {
            insFreq[j][0] = insFreqRef[j];
            delFreq[j][0] = delFreqRef[j];
        }

        for (int i = 0; i < orgIds.length; i++) {
            String[][] arrOrg = database.getOrganisms(new String[]{orgIds[i]});
            Organism org = new Organism(arrOrg[0][0], arrOrg[0][1]);
            ArrayList orgList = new ArrayList();
            orgList.add(org);

            String[][] orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(new String[]{org.getID()}, pathwayIds);

            MetabolicNetwork orgNetwork = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
                    org.getName(), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            orgNetwork.markCoreReactions(coreReactions);
            orgNetwork.addInactiveReactions(refNetwork.getReactions());
            orgNetwork.removeInactiveMetabolites();
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);

            orgNetwork.setupNetworkMatrices();
            orgNetwork.setupNetworkSequence();

            // Get core and prohibited edges
            ArrayList coreEdges = MetabolicNetwork.createDirectedReactions(coreReactions);
            ArrayList prohibEdges = MetabolicNetwork.createDirectedReactions(prohibReactions);
            lstNetworks = networkSimulator.simulateNetworkEvolution(orgNetwork, NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, nIter, coreEdges, prohibEdges, insRate, delRate, updateInterval, 1.0);
            // variables to hold summary data 
            int[] insFreqOrg = new int[orgNetwork.getDirectedReactions().size()]; // number of insertions for each edge
            int[] delFreqOrg = new int[orgNetwork.getDirectedReactions().size()]; // number of deletions for each edge

            NetworkEvolver.summariseEvolution(lstNetworks, insFreqOrg, delFreqOrg, null, null);

            for (int j = 0; j < insFreqOrg.length; j++) {
                insFreq[j][i + 1] = insFreqOrg[j];
                delFreq[j][i + 1] = delFreqOrg[j];
            }
        }

        System.out.println("Insertion Frequency:");
        System.out.print("ID" + "\t" + "REF" + "\t");
        for (int i = 0; i < orgIds.length; i++)
            System.out.print(orgIds[i] + "\t");
        System.out.println();
        for (int i = 0; i < insFreq.length; i++) {
            int[] rxnInsFreq = insFreq[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < rxnInsFreq.length; j++) {
                System.out.print(rxnInsFreq[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println("Deletion Frequency:");
        System.out.print("ID" + "\t" + "REF" + "\t");
        for (int i = 0; i < orgIds.length; i++)
            System.out.print(orgIds[i] + "\t");
        System.out.println();
        for (int i = 0; i < delFreq.length; i++) {
            int[] rxnDelFreq = delFreq[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < rxnDelFreq.length; j++) {
                System.out.print(rxnDelFreq[j]);
                System.out.print("\t");
            }
            System.out.println();
        }
    }

    static public void compareEvolutionForDifferentRates(int mode, String[] orgIds, String[] pathwayIds, double dependenceProbability) {

        int nIter = 1000;
        int updateInterval = 100;
        double[] insRates = new double[]{0.2, 0.4, 0.6, 0.8, 1.0};
        double[] delRates = new double[]{1.0, 0.8, 0.6, 0.4, 0.2};

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
        String ignoreFile = "ignoreCompoundsTCA.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));


        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        int[][] insFreq = new int[refNetwork.getDirectedReactions().size()][insRates.length]; // number of insertions for each edge
        int[][] delFreq = new int[refNetwork.getDirectedReactions().size()][insRates.length]; // number of deletions for each edge
        int[][] edgeCount = new int[nIter + 1][insRates.length]; // number of deletions for each edge

        MetabolicNetwork theNetwork;
        if (orgIds == null || orgIds.length == 0)
            theNetwork = refNetwork;
        else {
            String[][] arrOrg = database.getOrganisms(orgIds);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds, pathwayIds);

            theNetwork = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            theNetwork.addInactiveReactions(refNetwork.getReactions());
            theNetwork.removeInactiveMetabolites();
            theNetwork.getMetabolites().removeAll(ignoreMetabolites);

            if (mode != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                theNetwork.setupNetworkMatrices();

            theNetwork.setupNetworkSequence();
        }

        for (int i = 0; i < insRates.length; i++) {
            double insRate = insRates[i];
            double delRate = delRates[i];

            NetworkSimulator networkSimulator = new NetworkSimulator();
            ArrayList lstNetworks = networkSimulator.simulateNetworkEvolution(theNetwork, mode, nIter, new ArrayList(), new ArrayList(), insRate, delRate, updateInterval, dependenceProbability);

            // variables to hold summary data 
            int[] insFreqOrg = new int[theNetwork.getDirectedReactions().size()]; // number of insertions for each edge
            int[] delFreqOrg = new int[theNetwork.getDirectedReactions().size()]; // number of deletions for each edge
            int[] edgeCountOrg = new int[nIter + 1]; // number of deletions for each edge

            NetworkEvolver.summariseEvolution(lstNetworks, insFreqOrg, delFreqOrg, null, edgeCountOrg);

            for (int j = 0; j < insFreqOrg.length; j++) {
                insFreq[j][i] = insFreqOrg[j];
                delFreq[j][i] = delFreqOrg[j];
            }
            for (int j = 0; j < edgeCountOrg.length; j++) {
                edgeCount[j][i] = edgeCountOrg[j];
            }

        }

        System.out.println("Insertion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 0; i < insRates.length; i++) {
            System.out.print(insRates[i] / delRates[i]);
            System.out.print("\t");
        }
        System.out.println();
        for (int i = 0; i < insFreq.length; i++) {
            int[] rxnInsFreq = insFreq[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < rxnInsFreq.length; j++) {
                System.out.print(rxnInsFreq[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println("Deletion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 0; i < insRates.length; i++) {
            System.out.print(insRates[i] / delRates[i]);
            System.out.print("\t");
        }
        System.out.println();
        for (int i = 0; i < delFreq.length; i++) {
            int[] rxnDelFreq = delFreq[i];
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < rxnDelFreq.length; j++) {
                System.out.print(rxnDelFreq[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println("Edge Count:");
        System.out.print("Iter" + "\t");
        for (int i = 0; i < insRates.length; i++) {
            System.out.print(insRates[i] / delRates[i]);
            System.out.print("\t");
        }
        System.out.println();
        for (int i = 0; i < edgeCount.length; i++) {
            int[] nwEdgeCount = edgeCount[i];
            System.out.print(i);
            System.out.print("\t");
            for (int j = 0; j < nwEdgeCount.length; j++) {
                System.out.print(nwEdgeCount[j]);
                System.out.print("\t");
            }
            System.out.println();
        }

    }

    public void simulateKOCascading(MetabolicNetwork startNetwork, ArrayList coreEdges, ArrayList prohibEdges,
            int evolutionModel, double insRate, double delRate, ArrayList ignoreMetabolites, double dependenceProbability) {

        int neighbourCutoff = 1;
        MetabolicNetwork theNetwork = startNetwork.clone();
        theNetwork = theNetwork.filterNetworkByNeighbourhood(neighbourCutoff, ignoreMetabolites);
        theNetwork.setDirectedReactions(MetabolicNetwork.createDirectedReactions(theNetwork.getReactions()));
        theNetwork.setupNetworkSequence();
//        // find the differences in the two networks
//        ArrayList lstDifferences = MetabolicNetwork.findDifferences(theNetwork, startNetwork);
//        if (lstDifferences.size() > 0) {
//            // print the ids
//            System.out.print("Initial Filtering: ");
//            Iterator it = lstDifferences.iterator();
//            while (it.hasNext()) {
//                DirectedReaction reaction = (DirectedReaction) it.next();
//                System.out.print(reaction.getID() + "\t");
//            }
//            System.out.println();
//        }
        
        boolean first = true;
        while (true) {
            // select a random hyperedge for KO
            double[] edgeRates = new NetworkEvolver().calculateHyperedgeChangeRates(theNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);
            Byte[] rxnSeq = theNetwork.getReactionSequence();
            for (int i = 0; i < rxnSeq.length; i++) {
                // remove rates for insertable edges (only KOs required)
                if (rxnSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT))
                    edgeRates[i] = 0.0;
            }
            // normalise the rates to get probabilities
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

            // select a hyperedge for KO based on edge probabilities
            DirectedReaction hyperedge = (DirectedReaction) theNetwork.getDirectedReactions().get(Utilities.randsample(edgeProbabilities));
            // update the network
            MetabolicNetwork newNetwork = theNetwork.clone();
            newNetwork.deleteReaction(hyperedge.getReactionId());
            // apply neighbourhood filtering
            newNetwork = newNetwork.filterNetworkByNeighbourhood(neighbourCutoff, ignoreMetabolites);
            newNetwork.setDirectedReactions(MetabolicNetwork.createDirectedReactions(newNetwork.getReactions()));
            newNetwork.setupNetworkSequence();

            if (first) {
                System.out.println("Beginning KO");
                first = false;
            }
            // print the ids
            System.out.println("\tKO Edge: " + hyperedge.getID());
            System.out.print("\tCascade Effect: ");
            
            // find the differences in the two networks
            ArrayList lstDifferences = MetabolicNetwork.findDifferences(theNetwork, newNetwork);
            // stop if no additional change
            if (lstDifferences.size() == 1) {
                System.out.println("\t(None) ");
                System.out.println();
                break;
            } else {
                Iterator it = lstDifferences.iterator();
                while (it.hasNext()) {
                    DirectedReaction reaction = (DirectedReaction) it.next();
                    if (reaction.equals(hyperedge))
                        continue;
                    System.out.print(reaction.getID() + "\t");
                }
                System.out.println();
                System.out.println();
            }

            theNetwork = newNetwork;
        } // end while true
    }
}
