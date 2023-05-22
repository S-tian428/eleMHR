package eleMHRdemo;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class eleMHR {
    //For five elements in a cancer type and write out
    public static void main(String[] args) throws IOException {
        long startTimes = System.currentTimeMillis();

        //Read the mutation file of elements  1ms
        //Set current path
        String root = "path\\data\\test_file\\element_mutations\\";
        File dir = new File(root);
        //Gets all file names in the path
        String[] names = dir.list();
        for (int i = 0; i < names.length; i++) {
            String filename = names[i];
            //test
            //String filename = names[1];
            // System.out.println(filename);

            //BufferedWriter
            String WriterFile = "path\\data\\test_file\\Java_element_result\\";
            BufferedWriter bw = new BufferedWriter(new FileWriter(WriterFile + filename));
            //Write out title first
            bw.write("gene" + "\t" + "p" + "\t" + "nMut" + "\t" + "es.max" + "\t" + "maxpos" + "\t" + "es.min" + "\t" + "minpos" + "\t" + "es" + "\n");
            //Get start site for each gene    4ms
            LinkedHashMap<String, Integer> geneStart = getStarts(root, filename);
            //System.out.println(geneStart);
            //Get end site for each gene    1ms
            LinkedHashMap<String, Integer> geneEnd = getEnds(root, filename);
            //System.out.println(geneEnd);
            //Get the mutation locations in each gene  1ms
            LinkedHashMap<String, ArrayList<Integer>> genePositions = getPositions(root, filename);
            //System.out.println(genePositions);
            //For each gene in the genePositions map
            Set<String> genes = genePositions.keySet();
            for (String gene : genes) {
                // String gene = "CDH23";
                //Obtain the set of gene mutation location  0ms
                ArrayList<Integer> mutPositions = genePositions.get(gene);
                //Get the total number of mutations in gene 0ms
                //System.out.println(mutPositions);
                int Nm = mutPositions.size();
                // System.out.println(Nm);
                //gene start site 0ms
                int start = geneStart.get(gene);
                //System.out.println(start);
                //gene end site 0ms
                int end = geneEnd.get(gene);
                //System.out.println(end);
                //gene length 0ms
                int length = end - start + 1;
                // System.out.println(length);
                // Generate a list of consecutive integers starting with start and ending with end 49ms
                List<Integer> range = IntStream.rangeClosed(start, end).boxed().collect(Collectors.toList());
                // System.out.println(range);
                //es.true 63ms
                double[] resultTrue = getResult(mutPositions, length, range, start);
                //Iterate resultTrue
         /*for (int i = 0; i < resultTrue.length; i++) {
            if (i == resultTrue.length-1 ){
                System.out.println(r
                esultTrue[i]);
            }else {
                System.out.print(resultTrue[i] + " ");
            }
         }*/

                //Generating random mutations for Nm*10000 times
                int times = 1000;
                ArrayList<Integer> list = mutPosRandom(range, Nm, times);//6毫秒
                //Generate a list and put esRandom result into it to convenient for subsequent calculation of mean value and standard deviation
                //Convert the esRandoms collection to thread-safe
                List<Double> esRandoms = Collections.synchronizedList(new LinkedList<>());
                //// take Nm of them one time and put them into list    7ms
                List<List<Integer>> lists = TrackSubpackage(list, Nm);
                //System.out.println(list);
                //Carry out the getResult
                lists.parallelStream().
                        forEachOrdered(integers -> esRandoms.add(getResult(integers, length, range, start)[4]));
            /*slow speed
        for (List<Integer> integers : lists) {
            double[] esRandom = getResult(integers, length, range,start);
            esRandoms.add(esRandom[4]);
        }*/
                // System.out.println(esRandoms);
                //Calculate the mean and standard deviation of ES random set
                //mean,3ms
                double avg = esRandoms.stream().collect(Collectors.averagingDouble(Double::doubleValue));
                // System.out.println(avg);
                //standard deviation 0ms
                double standardDeviation = getStandardDeviation(esRandoms, avg);
                //System.out.println(standardDeviation);
                //Calculate  p values based on the mean and variance of es.random
                double trueES = resultTrue[4];
                //p 1ms
                double p = getNumber(trueES, esRandoms) / times;
                //  System.out.println(p);
                //Output the result 1ms
                String str = gene + "\t" + p + "\t" + Nm + "\t" + resultTrue[0] + "\t" +
                        resultTrue[1] + "\t" + resultTrue[2] + "\t" + resultTrue[3] + "\t" + resultTrue[4] + "\n";
                //System.out.println(str);
                bw.write(str);
                System.out.println(gene);
            }
            bw.close();
            long endTimes = System.currentTimeMillis();
            System.out.println(((endTimes - startTimes)));
            // System.out.println(i);
            // System.out.println(esRandoms);
        }


    }

    //Functions
    //Function1:Input BufferedReader object to get the Map of element start sites for each gene
    public static LinkedHashMap<String, Integer> getStarts(String root, String filename) throws IOException {
        FileReader fr = new FileReader(root + filename);
        BufferedReader br = new BufferedReader(fr);
        //Delete the first line
        br.readLine();
        //Create linkedHashMap to record unique gene-start(Integer) object
        LinkedHashMap<String, Integer> geneStart = new LinkedHashMap<>();
        String line = null;
        while ((line = br.readLine()) != null) {
            String[] arrStr = line.split("\\t");
            //gene name(Column 4)
            String gene = arrStr[3];
            //element start sites（Column 1）
            int start = Integer.parseInt(arrStr[0]);
            //Store gene and start site into Map
             geneStart.put(gene, start);
        }
        br.close();
        return geneStart;
    }

    //Function2:Input BufferedReader object to get the Map of element end sites for each gene
    public static LinkedHashMap<String, Integer> getEnds(String root, String filename) throws IOException {
        FileReader fr = new FileReader(root + filename);
        BufferedReader br = new BufferedReader(fr);
        //Delete the first line
        br.readLine();
        //Create linkedHashMap to record unique gene-end(Integer) object
        LinkedHashMap<String, Integer> getEnd = new LinkedHashMap<>();
        String line = null;
        while ((line = br.readLine()) != null) {
            String[] arrStr = line.split("\\t");
            //gene name(Column 4)
            String gene = arrStr[3];
            // System.out.println(gene);
            //element end sites（Column 2）
            int end = Integer.parseInt(arrStr[1]);
            //Store gene and end site into Map
            getEnd.put(gene, end);
        }
        br.close();
        return (getEnd);
    }

    //Function3:Input BufferedReader object to get the Map of mutation sites for each gene
    public static LinkedHashMap<String, ArrayList<Integer>> getPositions(String root, String filename) throws IOException {
        FileReader fr = new FileReader(root + filename);
        BufferedReader br = new BufferedReader(fr);
        //Delete the first line
        br.readLine();

        String line = null;
        //Create Map with unique gene-multiple mutation sites
        LinkedHashMap<String, ArrayList<Integer>> mutPositions = new LinkedHashMap<>();
        //To avoid null calls, first read a row and store the first gene and the first mutation position into the Map
        line = br.readLine();
        //The first line produces a array
        String[] arrStrPre = line.split("\\t");
        // Get the gene for the first row
        String gene1 = arrStrPre[3];
        //Get the mutation position for the first row and put it into the list
        ArrayList<Integer> position1 = new ArrayList<>();
        //Mutation locations（Column 6）
        position1.add(Integer.parseInt(arrStrPre[5]));
        //Store the entire key-value pair into the Map
        mutPositions.put(gene1, position1);
        //For each line
        while ((line = br.readLine()) != null) {
            String[] arrStr = line.split("\\t");
            //genes
            String gene = arrStr[3];
            // System.out.println(gene);
            //mutation position
            int position = Integer.parseInt(arrStr[5]);
            //Get all keys set
            Set<String> keys = mutPositions.keySet();
            //Determine if all keys in the current map collection contain the current gene
            if (keys.contains(gene)) {
                //If contains, add it
                mutPositions.get(gene).add(position);
            } else {
                //If not,create a new list,store the mutation site of this gene
                ArrayList<Integer> value = new ArrayList<>();
                value.add(position);
                //And put the key-value pairs of the new gene into the Map
                mutPositions.put(gene, value);
            }

        }
        br.close();
        return (mutPositions);
    }


    //Function4: Count the number of mutations occurred at each position，and put it into a Map
    //0ms
    public static HashMap<Integer, Double> mutPosTimes(List<Integer> mutPosR) {
        HashMap<Integer, Double> hm = new HashMap<>();
        //Iterate list
        //position:mutation location
        for (Integer position : mutPosR) {
            //Whether the current site exists in the map
            if (hm.containsKey(position)) {
                //Yes
                hm.put(position, hm.get(position) + 1);

            } else {
                //No
                hm.put(position, 1.0);
            }
        }
        return hm;
    }

    //Function5:(value of Map) * inc
    //Multiply the same value, the order of the keys doesn't matter，so use the HashMap
    //0ms
    public static HashMap<Integer, Double> incMap(HashMap<Integer, Double> mutPosT, double inc) {
        //Iterate over the HashMap, extract the value and multiply it by inc
        Set<Integer> keys = mutPosT.keySet();
        for (Integer key : keys) {
            Double value = mutPosT.get(key) * inc;
            mutPosT.put(key, value);
        }
        return mutPosT;
    }

    //Function6:Make the consecutive set of start to end positions as KEY,-dec as value,put them into map
    //The keys need to be ordered, so use LinkedHashMap
    public static LinkedHashMap<Integer, Double> decMap(List<Integer> range, double dec) {

        LinkedHashMap<Integer, Double> allPosDec = new LinkedHashMap<>();
        //slower:range.parallelStream.forEach(position ->  allPosDec.put(position, -dec));
        //slower：range.forEach(position ->  allPosDec.put(position, -dec));
        for (Integer position : range) {
            allPosDec.put(position, -dec);
        }
        return allPosDec;
    }

    //Function7:Calculate the cumulative sum of values in Map,
    //calculate the maximum value, minimum value and the difference between them in cumulative sum list
    public static double[] getES(LinkedHashMap<Integer, Double> seqIncDec, int start) {
        //Take all values from Map to a set
        List<Double> valuesList = new ArrayList<>();
        Set<Map.Entry<Integer, Double>> entries = seqIncDec.entrySet();
        Iterator<Map.Entry<Integer, Double>> iterator = entries.iterator();
        //Get the first value
        double temp = entries.iterator().next().getValue();
        while (iterator.hasNext()) {
            //Start from the second value
            double value = iterator.next().getValue();
            //Add value to the first one
            temp = value + temp;
            //assign temp to the value
            valuesList.add(temp);
        }
        //计算集合的最大值最小值以及最大值最小值对应的基因组位置
        //Calculate the maximum and minimum values and the genome positions for them
        //  System.out.println(valuesList);
        double max = Collections.max(valuesList);
        // System.out.println(max);
        double min = Collections.min(valuesList);
        //  System.out.println(min);
        //The genomic position for the maximum value
        double maxPos = valuesList.indexOf(max) + start;
        // System.out.println(maxPos);
        //The genomic position for the minimum  value
        double minPos = valuesList.indexOf(min) + start;
        // System.out.println(minPos);
        //Calculate the difference between the max and min
        double es = max - min;
        double[] trueResult = {max, maxPos, min, minPos, es};
        return trueResult;

    }

    //Function8:Calculate the true and random es
    public static double[] getResult(List<Integer> mutPos, int length, List<Integer> range, int start) {
        //inc(1/The sum of all the mutations occurred)
        double inc = 1.0 / mutPos.size();
        //System.out.println(inc);
        //Count the number of mutations at each location 0ms
        HashMap<Integer, Double> mutPosT = mutPosTimes(mutPos);
        //System.out.println(mutPosT);
        //Multiply the number of mutations at each location by inc 0ms
        HashMap<Integer, Double> mutPosTInc = incMap(mutPosT, inc);
        // System.out.println(mutPosT);
        //dec(1/(gene length - the number of locations where mutations occurred))
        double dec = 1.0 / (length - mutPosT.size());
        //System.out.println(dec);
        //生成一个有序的map,键为基因起始位置和终止位置的连续集合，
        //Generate an ordered map, whose keys are a continuous set of where the genes start and end
        // values is -dec     //20ms
        LinkedHashMap<Integer, Double> seqIncDec = decMap(range, dec);
        // System.out.println(seqIncDec);
        //Merge the two maps, replace -dec at the mutate location with the inc*mutations number 0ms
        seqIncDec.putAll(mutPosTInc);
        //System.out.println(seqIncDec);
        //Get the cumulative sum of values in the combined map,
        // and the max,min,max position,min position, and (max-min) of result list
        //dozens ms
        double[] es = getES(seqIncDec, start);

        //  System.out.println("es result：");
       /* for (int i = 0; i < es.length; i++) {
            if (i == es.length - 1) {
                System.out.println(es[i]);
            } else
                System.out.print(es[i] + " ");
        }*/
        return es;
    }

    // 9、Randomly generate Nm mutation sites according to the continuous set of genes sequence,and put them into Map
    public static ArrayList<Integer> mutPosRandom(List<Integer> range, int mutNum, int times) {

        ArrayList<Integer> mutPos = new ArrayList<>();
        //Random mutNum*times times,and put the random locations into the collection
        //seed：Same number of random seeds to make sure every randomize result is the same,in turn ensure the repeatability
        Random r = new Random(1);
        for (int i = 0; i < mutNum * times; i++) {
            //The random index of the List is produced, not the value
            //There is no need to ensure non-duplication
            int randomIndex = r.nextInt(range.size());
            mutPos.add(range.get(randomIndex));
        }
        return mutPos;
    }

    //Function10：Extract Nm sets at a time from the entire mutNum*times
    public static List<List<Integer>> TrackSubpackage(ArrayList<Integer> objects, int everySize) {
        List<List<Integer>> packageList = new ArrayList<>();
        //packages as a temporary set reader, and release it after filling mutNum at a time
        List<Integer> packages = new ArrayList<>();
        int size = objects.size();
        for (int i = 0; i < size; i++) {
            packages.add(objects.get(i));
            if (packages.size() == everySize || i == objects.size() - 1) {
                List<Integer> isp = new ArrayList<>();
                isp.addAll(packages);
                packageList.add(isp);
                packages.clear();
            }
        }
        return packageList;
    }

    //Function11：Calculate the standard deviation of the set
    //The denominator of the standard deviation formula is (n - 1)
    public static double getStandardDeviation(List<Double> esRandoms, double avg) {
        double dVar = 0;
        double size = esRandoms.size();
        for (int i = 0; i < size; i++) {
            //Calculate the variance
            dVar += (esRandoms.get(i) - avg) * (esRandoms.get(i) - avg);
        }

        return Math.sqrt(dVar / (size - 1));
    }

    //Function12：Count the number of the random es values that greater than the true
    public static double getNumber(double trueES, List<Double> esRandoms) {
        double number = 0;
        for (int i = 0; i < esRandoms.size(); i++) {
            if (esRandoms.get(i) >= trueES) {
                number++;
            }
        }
        return number;
    }


}

