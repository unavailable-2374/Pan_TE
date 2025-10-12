import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.nio.file.*;

public class ParallelKmersFeaturesCollector {
    
    private static final int NUM_THREADS = Math.max(1, Runtime.getRuntime().availableProcessors() * 4 / 5);
    private static final List<String> NUCLEOTIDES = Arrays.asList("A", "C", "T", "G");
    
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        
        File currentDir = new File(new File(".").getAbsolutePath());
        String currDirPath = currentDir.getCanonicalPath();
        
        String listFilePath = currDirPath + "/list.txt";
        List<String> sequenceIds = Files.readAllLines(Paths.get(listFilePath));
        
        String outputFeatureFilePath = currDirPath + "/feature_file.csv";
        
        System.out.println("Processing " + sequenceIds.size() + " sequences using " + NUM_THREADS + " threads");
        
        // Initialize feature map template
        Map<String, Integer> featureTemplate = initializeFeatureMap();
        List<String> orderedFeatures = getOrderedFeatureList();
        
        // Write header
        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFeatureFilePath))) {
            writer.println(String.join(",", orderedFeatures));
            
            // Process sequences in parallel batches
            ExecutorService executor = Executors.newFixedThreadPool(NUM_THREADS);
            List<Future<String>> futures = new ArrayList<>();
            AtomicInteger processedCount = new AtomicInteger(0);
            
            for (String sequenceId : sequenceIds) {
                Future<String> future = executor.submit(() -> {
                    try {
                        String featureLine = processSequence(sequenceId, currDirPath, featureTemplate, orderedFeatures);
                        int completed = processedCount.incrementAndGet();
                        if (completed % 100 == 0) {
                            System.out.println("Processed " + completed + "/" + sequenceIds.size() + " sequences");
                        }
                        return featureLine;
                    } catch (Exception e) {
                        System.err.println("Error processing " + sequenceId + ": " + e.getMessage());
                        return null;
                    }
                });
                futures.add(future);
            }
            
            // Collect results in order
            for (Future<String> future : futures) {
                String featureLine = future.get();
                if (featureLine != null) {
                    writer.println(featureLine);
                }
            }
            
            executor.shutdown();
        }
        
        System.out.println("Feature collection completed successfully!");
    }
    
    private static String processSequence(String sequenceId, String currDirPath, 
                                        Map<String, Integer> featureTemplate, 
                                        List<String> orderedFeatures) throws IOException {
        
        // Create a copy of the feature map for this sequence
        Map<String, Integer> featureMap = new HashMap<>(featureTemplate);
        
        // File paths for k-mer results
        String twoMerFile = currDirPath + "/kanalyze-2.0.0/output_data/2mer/" + sequenceId + ".txt";
        String threeMerFile = currDirPath + "/kanalyze-2.0.0/output_data/3mer/" + sequenceId + ".txt";
        String fourMerFile = currDirPath + "/kanalyze-2.0.0/output_data/4mer/" + sequenceId + ".txt";
        
        // Read k-mer counts from files
        readKmerFile(twoMerFile, featureMap);
        readKmerFile(threeMerFile, featureMap);
        readKmerFile(fourMerFile, featureMap);
        
        // Generate feature line
        StringBuilder features = new StringBuilder();
        for (int i = 0; i < orderedFeatures.size(); i++) {
            if (i > 0) features.append(",");
            features.append(featureMap.get(orderedFeatures.get(i)));
        }
        
        return features.toString();
    }
    
    private static void readKmerFile(String filePath, Map<String, Integer> featureMap) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] tokens = line.split("\t");
                if (tokens.length >= 2) {
                    String kmer = tokens[0];
                    Integer count = Integer.parseInt(tokens[1]);
                    featureMap.put(kmer, count);
                }
            }
        }
    }
    
    private static Map<String, Integer> initializeFeatureMap() {
        Map<String, Integer> featureMap = new HashMap<>();
        
        // Generate 2-mers
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                featureMap.put(n1 + n2, 0);
            }
        }
        
        // Generate 3-mers
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                for (String n3 : NUCLEOTIDES) {
                    featureMap.put(n1 + n2 + n3, 0);
                }
            }
        }
        
        // Generate 4-mers
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                for (String n3 : NUCLEOTIDES) {
                    for (String n4 : NUCLEOTIDES) {
                        featureMap.put(n1 + n2 + n3 + n4, 0);
                    }
                }
            }
        }
        
        return featureMap;
    }
    
    private static List<String> getOrderedFeatureList() {
        List<String> features = new ArrayList<>();
        
        // Add 2-mers in order
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                features.add(n1 + n2);
            }
        }
        
        // Add 3-mers in order
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                for (String n3 : NUCLEOTIDES) {
                    features.add(n1 + n2 + n3);
                }
            }
        }
        
        // Add 4-mers in order
        for (String n1 : NUCLEOTIDES) {
            for (String n2 : NUCLEOTIDES) {
                for (String n3 : NUCLEOTIDES) {
                    for (String n4 : NUCLEOTIDES) {
                        features.add(n1 + n2 + n3 + n4);
                    }
                }
            }
        }
        
        return features;
    }
}